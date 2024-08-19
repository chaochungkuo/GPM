import logging
import sys
from pathlib import Path
from typing import Optional, Union

import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from celltypist import models as ctypist_models
from pandas import DataFrame

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             Minor annotatio functions                                                                                    ###
### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###


def cell_typist_annotate(
    adata: AnnData, models: list[str], inplace=True, COUNTS_LAYER="counts"
) -> AnnData:
    if len(models) == 0:
        raise ValueError("The models list are empty, enter valid model names.")

    all_models = ctypist_models.models_description().model.to_list()

    for model in models:
        if model not in all_models:
            raise ValueError("{model} not found in supported cell typist models.")

    ctypist_models.download_models(force_update=True, model=models)

    adata_celltypist = adata.copy()
    adata_celltypist.X = adata.layers[COUNTS_LAYER]
    sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)
    sc.pp.log1p(adata_celltypist)
    adata_celltypist.X = adata_celltypist.X.toarray()

    for model in models:
        loaded_model = ctypist_models.Model.load(model=model)
        predictions = celltypist.annotate(
            adata_celltypist, model=loaded_model, majority_voting=True
        )
        predictions_adata = predictions.to_adata()
        adata.obs["celltypist_" + model + "_annotation"] = predictions_adata.obs.loc[
            adata.obs.index, "majority_voting"
        ]
        adata.obs["celltypist_" + model + "_conf_score"] = predictions_adata.obs.loc[
            adata.obs.index, "conf_score"
        ]
    if not inplace:
        return adata


def get_marker_genes(marker_gene_path: Path, adata: AnnData) -> pd.DataFrame:
    """
    Get marker genes for each cluster
    """
    df = pd.read_csv("./marker_genes.csv", header=None)
    columns: list[str] = ["cell_type", "gene"]
    for i in range(2, df.shape[1]):
        columns.append(f"col{i}")
    df.columns = columns
    df: pd.DataFrame = df[df.gene.isin(adata.var_names)]

    if len(df) == 0:
        raise ValueError("No marker genes found in the dataset")


def score_markers(markers_df, adata: AnnData, annotation_threshold: float) -> None:
    gf = markers_df.groupby("cell_type")
    for type in markers_df.cell_type.unique():
        gene_list = gf.get_group(type).gene.to_list()
        print(gene_list)
        sc.tl.score_genes(adata, gene_list, score_name=f"{type}_score")

    adata.obs["marker_annotation"] = "unassigned"

    for cell_type in markers_df.cell_type.unique():
        adata.obs["marker_annotation"] = np.where(
            adata.obs[f"{cell_type}_score"] > annotation_threshold,
            cell_type,
            adata.obs["marker_annotation"],
        )

    # Mark cells where multiple idents are above the threshold as 'ambiguous'
    cols_to_susbet: list[str] = [
        f"{cell_type}_score" for cell_type in markers_df.cell_type.unique()
    ]
    adata.obs["marker_annotation"] = np.where(
        (adata.obs[cols_to_susbet] > annotation_threshold).sum(axis=1) > 1,
        "ambiguous",
        adata.obs["marker_annotation"],
    )


### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                             Majority voting of predictions                                                                               ###
### -------------------------------------------------------------------------------------------------------------------------------------------------------------------------###


# Over-clustering and majority voting functions are adapted from the cell-typist package by the Teichmann lab
# https://github.com/Teichlab/celltypist/blob/main/celltypist/classifier.py
def over_cluster(
    adata, resolution: Optional[float] = None, use_GPU: bool = False
) -> pd.Series:
    if use_GPU and "rapids_singlecell" not in sys.modules:
        logger.warn(
            "⚠️ Warning: rapids_singlecell is not installed but required for GPU running, will switch back to CPU"
        )
        use_GPU = False
    if "connectivities" not in adata.obsp:
        logger.info(
            "👀 Can not detect a neighborhood graph, will construct one before the over-clustering"
        )
        adata = sc.pp.pca(adata)
        adata = sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    else:
        logger.info(
            "👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it"
        )
    if resolution is None:
        if adata.n_obs < 5000:
            resolution = 5
        elif adata.n_obs < 20000:
            resolution = 10
        elif adata.n_obs < 40000:
            resolution = 15
        elif adata.n_obs < 100000:
            resolution = 20
        elif adata.n_obs < 200000:
            resolution = 25
        else:
            resolution = 30
    logger.info(f"⛓️ Over-clustering input data with resolution set to {resolution}")
    if use_GPU:
        sc.tl.leiden(adata, resolution=resolution, key_added="over_clustering")
    else:
        sc.tl.leiden(adata, resolution=resolution, key_added="over_clustering")
    return adata.obs.pop("over_clustering")


def majority_vote(
    predictions: pd.DataFrame,
    column: str,
    over_clustering: Union[list, tuple, np.ndarray, pd.Series, pd.Index],
    min_prop: float = 0,
) -> pd.DataFrame:
    if isinstance(over_clustering, (list, tuple)):
        over_clustering = np.array(over_clustering)
    logger.info("🗳️ Majority voting the predictions")
    votes = pd.crosstab(predictions[column], over_clustering)

    majority = votes.idxmax(axis=0).astype(str)
    freqs = (votes / votes.sum(axis=0).values).max(axis=0)
    majority[freqs < min_prop] = "Heterogeneous"
    majority: DataFrame = majority[over_clustering].reset_index()
    majority.index = predictions.index
    majority.columns = ["over_clustering", "majority_voting"]
    majority["majority_voting"] = majority["majority_voting"].astype("category")
    # predictions.predicted_labels = predictions.predicted_labels.join(majority)
    logger.info("✅ Majority voting done!")
    return majority
