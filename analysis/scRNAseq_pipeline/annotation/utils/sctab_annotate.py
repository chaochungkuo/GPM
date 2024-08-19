import argparse
from collections import OrderedDict
from os import path
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tomlkit
import torch
import yaml
from anndata import AnnData
from huggingface_hub import hf_hub_download
from scipy.sparse import csc_matrix
from sctab_utils import TabNet, dataloader_factory, streamline_count_matrix
from tqdm import tqdm
from util_funcs import majority_vote, over_cluster

###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             command line argument parsing                                                                               ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###

parser = argparse.ArgumentParser(
    prog="scTAB_annotate",
    description="Annotates an input scanpy object with scTAB model and saves the results to scTAB_annotation column",
    epilog="",
)


parser.add_argument(
    "-i",
    "--input",
    dest="filename",
    required=True,
    type=str,
    help="input path of scanpy file on disk.",
)
parser.add_argument(
    "-c",
    "--column",
    dest="column_name",
    default="scTAB_annotation",
    type=str,
    help="column to save the results to (default: scTAB_annotation)",
)

parser.add_argument(
    "-o",
    "--output",
    type=str,
    dest="ouput_path",
    default="",
    help="output path of the annotated scanpy file",
)

parser.add_argument(
    "-b",
    "--batchsize",
    dest="batchsize",
    default=2048,
    type=int,
    help="batch size for the data loader (default: 2048)",
)

parser.add_argument(
    "--config", dest="config", required=True, type=str, help="path to the config file"
)


args: argparse.Namespace = parser.parse_args()

if len(args.ouput_path) == 0:
    output_path: Path = args.filename
else:
    output_path: Path = Path(args.ouput_path)


###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             Config loading                                                                                              ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###

with open(args.config, "r") as f:
    config = tomlkit.parse(f.read())

ANALYSIS_DIR: str = config["basic"]["ANALYSIS_DIR"]
RESOURCE_PATH: str = config["basic"]["RESOURCE_PATH"]
DIR_SAVE: str = config["basic"]["DIR_SAVE"]
COUNTS_LAYER = config["normalization"]["COUNTS_LAYER"]

# print(
#     f"Resource path: {RESOURCE_PATH}",
#     f"Save path: {DIR_SAVE}",
#     f"filename: {args.filename}",
#     f"output_path {output_path}",
# )

###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             scTAB model functions                                                                                       ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###


def get_scTAB_resources() -> tuple[Path, Path, Path, Path]:
    weights: Path = hf_hub_download(
        "MohamedMabrouk/scTab",
        "val_f1_macro_epoch=41_val_f1_macro=0.847.ckpt",
        local_dir=path.join(RESOURCE_PATH, "scTAB"),
    )

    genes: Path = hf_hub_download(
        "MohamedMabrouk/scTab",
        "var.parquet",
        subfolder="merlin_cxg_2023_05_15_sf-log1p_minimal",
        local_dir=path.join(RESOURCE_PATH, "scTAB"),
    )

    hyperparams: Path = hf_hub_download(
        "MohamedMabrouk/scTab",
        "hparams.yaml",
        local_dir=path.join(RESOURCE_PATH, "scTAB"),
    )

    cell_type: Path = hf_hub_download(
        "MohamedMabrouk/scTab",
        "cell_type.parquet",
        subfolder="merlin_cxg_2023_05_15_sf-log1p_minimal/categorical_lookup",
        local_dir=path.join(RESOURCE_PATH, "scTAB"),
    )

    return weights, hyperparams, genes, cell_type


def scTAB_data_loader(adata: AnnData, genes_path: str, batchsize: int = 2048):
    genes_from_model: pd.DataFrame = pd.read_parquet(genes_path)

    # subset gene space only to genes used by the model
    adata = adata[:, adata.var.index.isin(genes_from_model.feature_name)]
    # pass the count matrix in csc_matrix to make column slicing efficient
    x_streamlined = streamline_count_matrix(
        csc_matrix(adata.X),
        adata.var.index,  # change this if gene names are stored in different column
        genes_from_model.feature_name,
    )
    loader = dataloader_factory(x_streamlined, batch_size=batchsize)

    return loader


def get_scTAB_model(weights_path: str, hyperparams_path: str) -> TabNet:
    # load checkpoint
    if torch.cuda.is_available():
        ckpt = torch.load(weights_path)
    else:
        # map to cpu if there is not gpu available
        ckpt = torch.load(
            weights_path,
            map_location=torch.device("cpu"),
        )

    tabnet_weights = OrderedDict()
    for name, weight in ckpt["state_dict"].items():
        if "classifier." in name:
            tabnet_weights[name.replace("classifier.", "")] = weight

    with open(hyperparams_path) as f:
        model_params = yaml.full_load(f.read())

    # initialzie model with hparams from hparams.yaml file
    tabnet = TabNet(
        input_dim=model_params["gene_dim"],
        output_dim=model_params["type_dim"],
        n_d=model_params["n_d"],
        n_a=model_params["n_a"],
        n_steps=model_params["n_steps"],
        gamma=model_params["gamma"],
        n_independent=model_params["n_independent"],
        n_shared=model_params["n_shared"],
        epsilon=model_params["epsilon"],
        virtual_batch_size=model_params["virtual_batch_size"],
        momentum=model_params["momentum"],
        mask_type=model_params["mask_type"],
    )

    # load trained weights
    tabnet.load_state_dict(tabnet_weights)
    # set model to inference mode
    tabnet.eval()

    return tabnet


def sf_log1p_norm(x):
    """Normalize each cell to have 10000 counts and apply log(x+1) transform."""

    counts = torch.sum(x, dim=1, keepdim=True)
    # avoid zero division error
    counts += counts == 0.0
    scaling_factor = 10000.0 / counts

    return torch.log1p(scaling_factor * x)


def scTAB_annotate() -> None:
    adata: AnnData = sc.read_h5ad(args.filename)
    validate_count_layer(adata, COUNTS_LAYER)
    weights, hyperparams, genes, cell_type = get_scTAB_resources()
    tabnet: TabNet[yaml.Any, yaml.Any] = get_scTAB_model(weights, hyperparams)
    loader = scTAB_data_loader(adata, genes, batchsize=2048)
    preds: list = []
    with torch.no_grad():
        for batch in tqdm(loader):
            # normalize data
            x_input = sf_log1p_norm(batch[0]["X"])
            logits, _ = tabnet(x_input)
            preds.append(torch.argmax(logits, dim=1).numpy())

    preds = np.hstack(preds)

    cell_type_mapping: pd.DataFrame = pd.read_parquet(cell_type)
    preds = cell_type_mapping.loc[preds]["label"].to_numpy()
    adata.obs[args.column_name] = pd.Categorical(preds)

    over_cluster_result: pd.Series = over_cluster(adata, use_GPU=False)
    majority_voting: pd.DataFrame = majority_vote(
        adata.obs, column=args.column_name, over_clustering=over_cluster_result
    )

    adata.obs[args.column_name + "_majority_voting"] = majority_voting[
        "majority_voting"
    ]

    # convert the columns to categorical
    adata.obs[args.column_name + "_majority_voting"] = pd.Categorical(
        adata.obs[args.column_name + "_majority_voting"]
    )
    adata.obs[args.column_name] = pd.Categorical(adata.obs[args.column_name])

    # Write the result to disk
    adata.write_h5ad(output_path)


def validate_count_layer(adata: AnnData, counts_layer: str = COUNTS_LAYER) -> None:
    if counts_layer == "X":
        adata.layers["counts"] = adata.X.copy()
        counts_layer = "counts"
    elif counts_layer in adata.layers.keys():
        adata.X = adata.layers[counts_layer].copy()
    else:
        raise ValueError("{counts_layer} layer can't be found in the object")


if __name__ == "__main__":
    scTAB_annotate()
