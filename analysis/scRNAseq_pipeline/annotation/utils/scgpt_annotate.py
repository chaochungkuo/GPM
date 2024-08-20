import argparse

# os.environ["CUDA_VISIBLE_DEVICES"] = "1"
# os.environ["WORLD_SIZE"] = "1"
import logging
import os
from os import path
from pathlib import Path

import faiss
import numpy as np
import scanpy as sc
import scgpt as scg
import tomlkit
from build_atlas_index_faiss import load_index, vote
from huggingface_hub import snapshot_download
from torch.cuda import is_available
from tqdm import tqdm

logging.basicConfig(level=logging.ERROR)

### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             command line argument parsing                                                                                ###
### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###


parser = argparse.ArgumentParser(
    prog="scGPT_annotate",
    description="Annotates an input scanpy object with scGPT CT model and saves the results to scGPT columns",
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
    default="scGPT_annotation",
    type=str,
    help="column to save the results to (default: scGPT_annotation)",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    dest="ouput_path",
    default="",
    help="output path for the annotated scanpy object",
)


parser.add_argument(
    "-b",
    "--batch-size",
    dest="batch_size",
    default=1024,
    type=int,
    help="batch size for embedding",
)


parser.add_argument(
    "--config", dest="config", required=True, type=str, help="path to the config file"
)


args = parser.parse_args()

if len(args.ouput_path) == 0:
    output_path = args.filename
else:
    output_path = args.ouput_path


with open(args.config, "r") as f:
    config = tomlkit.parse(f.read())

ANALYSIS_DIR: str = config["basic"]["ANALYSIS_DIR"]
RESOURCE_PATH: str = config["basic"]["RESOURCE_PATH"]
DIR_SAVE: str = config["basic"]["DIR_SAVE"]
COUNTS_LAYER = config["normalization"]["COUNTS_LAYER"]


def get_scGPT_resources() -> Path:
    model: Path = snapshot_download(
        "MohamedMabrouk/scGPT",
        local_dir=path.join(RESOURCE_PATH, "scGPT"),
    )
    return model


model_dir = get_scGPT_resources()
adata = sc.read_h5ad(args.filename)
# cell_type_key = "scGPT_metainfo"
gene_col = "index"
# adata.obs[cell_type_key] = 0

print(is_available())

###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                                         Gene Embeding                                                                                   ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
adata_embed = scg.tasks.embed_data(
    adata,
    model_dir,
    gene_col=gene_col,
    # obs_to_save=cell_type_key,  # optional arg, only for saving metainfo
    batch_size=args.batch_size,
    return_new_adata=True,
)
adata_embed = adata_embed.X

###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                                 Loading Index                                                                                           ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###

use_gpu = faiss.get_num_gpus() > 0
index, meta_labels = load_index(
    index_dir=path.join(model_dir, "cxg_faiss_index"),
    use_config_file=False,
    use_gpu=True,
)
print(f"Loaded index with {index.ntotal} cells")


###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            Annotating cells                                                                                             ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###


k = 50
# test with the first 100 cells
distances, idx = index.search(adata_embed, k)

predict_labels = meta_labels[idx]

voting = []
for preds in tqdm(predict_labels):
    voting.append(vote(preds, return_prob=False)[0])
voting = np.array(voting)
adata.obs[args.column_name] = voting

# convert the column to categorical
adata.obs[args.column_name] = adata.obs[args.column_name].astype("category")


###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             Writing results to disk                                                                                     ###
###-------------------------------------------------------------------------------------------------------------------------------------------------------------------------###

adata.write_h5ad(filename=output_path)
