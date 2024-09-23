import urllib.request
from collections.abc import Callable
from numbers import Number
from os import PathLike, path, system
from typing import Any

import numpy as np
import pandas as pd
import requests
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from scipy.sparse import issparse
from scipy.stats import median_abs_deviation

###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            Utility Functions                                                                               ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


# From: https://www.oreilly.com/library/view/python-cookbook/0596001673/ch04s16.html
def splitall(p) -> list[str]:
    """Splits a path into all its parts."""
    allparts = []
    while True:
        parts = path.split(p)
        if parts[0] == p:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == p:  # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            p = parts[0]
            allparts.insert(0, parts[1])
    return allparts


###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            Reading Functions                                                                               ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


def read_parsebio(data_path: PathLike) -> AnnData:
    """Reads the output of spipe parsebio pipeline and returns an AnnData object.

    Args:
        data_path (str): path to the parsebio output.

    Returns:
        AnnData: Return AnnData object containing .X, .obs, and .var components
    """

    adata: AnnData = sc.read_mtx(path.join(data_path, "count_matrix.mtx"))

    # reading in gene and cell data
    gene_data = pd.read_csv(path.join(data_path, "all_genes.csv"))
    cell_meta: pd.DataFrame = pd.read_csv(path.join(data_path, "cell_metadata.csv"))


    # find genes with nan values and filter
    gene_data: pd.DataFrame = gene_data[gene_data.gene_name.notnull()]
    not_na: pd.Index = gene_data.index
    not_na_ls: list[str] = not_na.to_list()

    # remove genes with nan values and assign gene names
    adata: AnnData = adata[:, not_na_ls]
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


###------------------------------------------------------------------------------------------------------------------------------------------------------------###
###                                                            QC Functions                                                                                    ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


def human2mouse(genes: list[str]) -> list[str]:
    """Converts human gene names to mouse gene names using the gprofiler API

    Args:
        genes (list[str]): a list of human gene names

    Returns:
        list[str]: converted list of mouse gene names. Failed conversions are replaced with pd.NA.
    """

    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler/api/orth/orth/",
        json={
            "organism": "hsapiens",
            "target": "mmusculus",
            "query": genes,
        },
        timeout=60,
    )
    df = pd.DataFrame(
        r.json()["result"],
    )
    return df.name.replace("N/A", pd.NA).dropna().to_list()

###---------------------------------------------------------------------------------------------------------------------------------###
###                                                            QC Functions                                                                                    ###
###------------------------------------------------------------------------------------------------------------------------------------------------------------###


def _compute_outliers(
    series: pd.Series,
    value: list | int | float,
    max_only: bool = False,
    log_transform: bool = False,
) -> pd.Series:
    """computes outliers for the given variable in the dataframe.

    Args:
        adata (AnnData): Input AnnData object.
        value (list | Number): value to use for outlier detection, if a list is provided, it is used as the lower and upper bound
        max_only (bool, optional): If True, only the upper bound is used for outlier detection. Defaults to False.
        log_transform (bool, optional): If True, the variable is log transformed before outlier detection. Defaults to False.
    Returns:
        df: the input dataframe with an additional column for the outliers
    """

    # Validate the input
    if not isinstance(value, (list, int, float)):
        raise ValueError(
            "Please provide a positive number of nmads or a list of length 2 for the lower and upper bound of the QC-variable."
        )

    if isinstance(value, list):
        min_val: float = value[0]
        max_val: float = value[1]

    if isinstance(value, (int, float)):
        if not value > 0:
            raise ValueError("Please provide a positive number of nmads.")

        if log_transform:
            series: np.ndarray = np.log1p(series)

        min_val: float = np.median(series) - (median_abs_deviation(series) * value)
        max_val: float = np.median(series) + (median_abs_deviation(series) * value)

    if max_only:
        outliers: pd.Series[bool] = series.gt(max_val)
    else:
        outliers: pd.Series[bool] = series.lt(min_val) | series.gt(max_val)
    outliers.fillna(False)

    return outliers

def compute_outliers(
    df: pd.DataFrame,
    qc_dict: dict[str, list | Number],
    max_only: list[str],
    log_transform: list[str],
) -> pd.DataFrame:

    missing_keys = [key for key in qc_dict.keys() if key not in df.columns]

    if len(missing_keys) > 0:
        raise KeyError(
            f"the following QC variables {','.join(missing_keys)} does not exist in the data, check the variable names again."
        )

    for key in qc_dict.keys():
        if f"{key}_outlier" not in df.columns:
            df[f"{key}_outlier"] = False

        if key in max_only:
            max_flag = True
        else:
            max_flag = False

        if key in log_transform:
            log_flag = True
        else:
            log_flag = False

        df[f"{key}_outlier"] = _compute_outliers(
            df[key], qc_dict[key], max_flag, log_flag
        )

    return df


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


def validate_qc_dict(dc: dict, df: pd.DataFrame) -> bool:
    if all(map(lambda x: isinstance(x, (list, Number)), dc.values())):
        return True
    elif all(map(lambda x: isinstance(x, (dict)), dc.values())):
        if df["sample"].unique() == dc.keys():
            return True

        raise ValueError(
            "The QC dictionary keys are either not covering all samples or the sample names are not matching the sample names in the dataframe."
        )
    else:
        return False


def GenomeInfoDB_fix(tmpdirname) -> None:
    # Workaround failure to install GenomeInfoDbData using pixi
    dn_path = path.join(tmpdirname, "GenomeInfoDbData_1.2.11.tar.gz")
    dn_url = "https://bioconductor.org/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
    urllib.request.urlretrieve(dn_url, filename=dn_path)
    system(f"R CMD INSTALL {dn_path}")


def create_panel_fig(
    *,
    total_plots,
    ncols=2,
    figsize=3,
    wspace=0.5,
    hspace=0.5,
) -> tuple[Figure, Any]:
    ncols = 2
    nrows = total_plots // ncols + total_plots % ncols

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(
            ncols * figsize + figsize * wspace * (ncols - 1),
            nrows * figsize + hspace * (nrows - 1),
        ),
    )
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    axes = axes.flatten()

    if len(axes) > total_plots:
        for i in range(total_plots, len(axes)):
            fig.delaxes(axes[i])
    return fig, axes
