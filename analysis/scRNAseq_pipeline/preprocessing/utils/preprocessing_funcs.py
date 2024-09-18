from numbers import Number
from os import path
from typing import Callable, Dict, List

import numpy as np
import pandas as pd
import requests
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse
from scipy.stats import median_abs_deviation

###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            Reading Functions                                                                               ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


# Configs
## Utility Functions
def get_sample_name(file_path: str, black_list: list[str], n=3):
    """ "Function to return probable sample name from a path, it recurselvy goes through the path and returns the first element not in the black list."""
    if n == 0:
        return ""

    tmp = path.basename(file_path)
    _d = path.dirname(file_path)

    if all(entry not in tmp for entry in black_list):
        return tmp
    else:
        res = get_sample_name(_d, black_list, n - 1)
    return res


def is_outlier(adata: AnnData, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def read_parsebio(data_path: str) -> AnnData:
    """Reads ParseBio

    Args:
        data_path (str): path to the parsebio output.

    Returns:
        AnnData: Return AnnData object containing .X, .obs, and .var components
    """

    adata = sc.read_mtx(data_path + "count_matrix.mtx")

    # reading in gene and cell data
    gene_data = pd.read_csv(data_path + "all_genes.csv")
    cell_meta = pd.read_csv(data_path + "cell_metadata.csv")

    # find genes with nan values and filter
    gene_data = gene_data[gene_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()

    # remove genes with nan values and assign gene names
    adata = adata[:, notNa]
    adata.var = gene_data
    adata.var.set_index("gene_name", inplace=True)
    adata.var.index.name = None
    adata.var_names_make_unique()

    # add cell meta data to anndata object
    adata.obs = cell_meta
    adata.obs.set_index("bc_wells", inplace=True)
    adata.obs.index.name = None
    adata.obs_names_make_unique()

    return adata


## Technology components

inputs: Dict[str, List | Callable] = {
    "10x": {
        "files": ["features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"],
        "black_list": ["filtered_feature_bc", "raw_feature_bc", "count", "outs"],
        "raw_name": "raw_feature_bc_matrix",
        "function": sc.read_10x_mtx,
    },
    "10_h5": {
        "files": ["*.h5"],
        "black_list": [],
        "raw_name": "",
        "function": sc.read_10x_h5,
    },
    "ParseBio": {
        "files": ["all_genes.csv", "cell_metadata.csv", "count_matrix.mtx"],
        "black_list": ["DGE_filtered", "DGE_unfiltered"],
        "function": read_parsebio,
        "raw_name": "",
    },
    "Singleron": {
        "files": ["features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"],
        "black_list": ["starsolo", "raw"],
        "function": sc.read_10x_mtx,
        "raw_name": "",
    },
}


qc_features_rules: Dict[str, List[str]] = {
    "human": {"mito": ["MT-"], "ribo": ["RBS", "RPL"], "hb": ["^HB[^(P)]"]},
    "mouse": {
        "mito": ["mt"],
        "ribo": ["Rps", "Rpl"],
        "hb": ["^Hb[^(p)]"],  # Validate this later
    },
}


###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            QC Functions                                                                                    ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


def human2mouse(genes: List[str]) -> List[str]:
    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler/api/orth/orth/",
        json={
            "organism": "hsapiens",
            "target": "mmusculus",
            "query": genes,
        },
    )
    df = pd.DataFrame(
        r.json()["result"],
    )
    return df.name.replace("N/A", pd.NA).dropna().to_list()


###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            Reading Functions                                                                               ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


## Technology components

inputs: Dict[str, List | Callable] = {
    "10x": {
        "files": ["features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"],
        "black_list": ["filtered_feature_bc", "raw_feature_bc", "count", "outs"],
        "raw_name": "raw_feature_bc_matrix",
        "function": sc.read_10x_mtx,
    },
    "10x_h5": {
        "files": ["*.h5"],
        "black_list": [],
        "raw_name": "",
        "function": sc.read_10x_h5,
    },
    "ParseBio": {
        "files": ["all_genes.csv", "cell_metadata.csv", "count_matrix.mtx"],
        "black_list": ["DGE_filtered", "DGE_unfiltered"],
        "function": read_parsebio,
    },
    "Singleron": {
        "files": ["features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"],
        "black_list": ["starsolo", "raw", "filtered"],
        "raw_name": "raw",
        "function": sc.read_10x_mtx,
    },
}


qc_features_rules: Dict[str, List[str]] = {
    "human": {"mito": ["MT-"], "ribo": ["RBS", "RPL"], "hb": ["^HB[^(P)]"]},
    "mouse": {
        "mito": ["mt"],
        "ribo": ["Rps", "Rpl"],
        "hb": ["^Hb[^(p)]"],  # Validate this later
    },
}


###------------------------------------------------------------------------------------------------------------------------------------###
###                                                            QC Functions                                                                                    ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


def _compute_outliers(
    df: pd.DataFrame,
    variable: str,
    value: List | Number,
    max_only: bool = False,
    log_transform: bool = False,
) -> pd.DataFrame:
    """computes outliers for the given variable in the dataframe

    Args:
        adata (AnnData): Input AnnData object
        variable: name of the column to use to for outlier detection
        value: value to use for outlier detection, if a list is provided, it is used as the lower and upper bound
    Returns:
        df: the input dataframe with an additional column for the outliers
    """

    if isinstance(value, list):

        min_val = value[0]
        max_val = value[1]
        max_only = False  # Make sure not to override the user input

    if isinstance(value, Number):
        if not value > 0:
            raise ValueError("Please provide a positive number of nmads")

        if log_transform:
            series = np.log1p(df[variable])
        else:
            series = df[variable]

        min_val = np.median(series) - (median_abs_deviation(series) * value)
        max_val = np.median(series) + (median_abs_deviation(series) * value)

    if max_only:
        df[f"{variable}_outlier"] = df[variable].gt(max_val)
    else:
        df[f"{variable}_outlier"] = df[variable].lt(min_val) | df[variable].gt(max_val)

    df[f"{variable}_outlier"].fillna(False)

    return df


def compute_outliers_var_bulk(
    df: pd.DataFrame,
    variable: str,
    value: List | Number,
    max_only: bool = False,
    log_transform: bool = False,
):
    
    # Validate the input
    if variable not in df:
        raise KeyError(
            "the provided QC variable does not exist in the data, check the variable names again."
        )

    if not isinstance(value, (list, Number)):
        raise ValueError(
            "Please provide a positive number of nmads or a list of length 2 for the lower and upper bound of the QC-variable."
        )

    if isinstance(value, Number):
        if not value > 0:
            raise ValueError("Please provide a positive number of nmads")


    return _compute_outliers(df, variable, value, max_only, log_transform) 





def compute_outlier_sample(
    adata: AnnData,
    variable: str,
    value: List | Number,
    sample: str,
    max_only: bool = False,
    log_transform: bool = False,
) -> None:
    adata.obs["outlier"] = False
    for sample in variables.keys():
        _compute_outlier_sample(adata, variables, sample)


def _compute_outlier_sample(adata: AnnData, variables: Dict[str, List], sample) -> None:
    sample_dict = variables[sample]
    sample_slice = adata.obs.query(f"sample == '{sample}'")

    for key in sample_dict.keys():
        if key not in sample_slice.columns:
            raise KeyError(
                "the provided QC variable does not exist in the data, check the variable names again."
            )

        if f"{key}_outlier" not in sample_slice.columns:
            sample_slice[f"{key}_outlier"] = False

        if isinstance(sample_dict[key], list):
            if len(sample_dict[key]) != 2:
                raise ValueError(
                    "Provide a list of length 2 for the lower and upper bound of the QC-variable."
                )

            min_val = sample_dict[key][0]
            max_val = sample_dict[key][1]

        if isinstance(sample_dict[key], Number):
            if not sample_dict[key] > 0:
                raise ValueError("Please provide a positive number of nmads")
            min_val = np.median(sample_slice[key]) - (
                median_abs_deviation(sample_slice[key]) * sample_dict[key]
            )
            max_val = np.median(adata.obs[key]) + (
                median_abs_deviation(sample_slice[key]) * sample_dict[key]
            )

        final_slice = sample_slice[key].lt(min_val) | sample_slice[key].gt(max_val)

        adata.obs.loc[adata.obs["sample"] == sample, f"{key}_outlier"] = final_slice
        adata.obs.loc[adata.obs["sample"] == sample, f"{key}_outlier"].fillna(False)

    # Aggregating outliers in the 'outlier' column
    adata.obs.loc[adata.obs["sample"] == sample, "outlier"] = adata.obs.loc[
        adata.obs["sample"] == sample, [f"{x}_outlier" for x in sample_dict.keys()]
    ].any(axis=1)


def get_keys(qc_dict):
    keys_list = []
    if len(qc_dict) > 0 and all(
        map(lambda x: isinstance(x, (list, Number)), qc_dict.values())
    ):
        return list(qc_dict.keys())

    if len(qc_dict) > 0 and all(
        map(lambda x: isinstance(x, (dict, Number)), qc_dict.values())
    ):
        for key in qc_dict.keys():
            keys_list = keys_list + list(qc_dict[key].keys())
        return list(set(keys_list))

    return []


###-----------------------------------------------------------------------------------------------------###
###                                       Normalization Functions                                       ###
###-----------------------------------------------------------------------------------------------------###


def get_var_features_num(adata: AnnData, variable_features: int | float) -> int:
    detected_gene_nu = len(adata.var_names)
    if variable_features <= 1:
        return int(detected_gene_nu * variable_features)
    else:
        return min(detected_gene_nu, variable_features)


def is_raw_counts(matrix) -> bool:

    if issparse(matrix):
        return matrix.count_nonzero() == matrix.astype("uint32").count_nonzero()
    else:
        return np.count_nonzero(matrix) == np.count_nonzero(matrix.astype("uint32"))


###-----------------------------------------------------------------------------------------------------###
###                                       Utility Functions                                       ###
###-----------------------------------------------------------------------------------------------------###


def select_n_uniform(length, n):
    if length < n:
        return list(range(length))
    step = length / n
    indices = [round(i * step) for i in range(n)]
    return indices
