from os import walk, path, mkdir, listdir
from scipy.stats import median_abs_deviation
import pandas as pd
import numpy as np
from anndata import AnnData
from typing import List, Dict, Callable
import scanpy as sc
import requests
from functools import reduce


# Configs
## Utility Functions
def get_sample_name(file_path: str, black_list: list[str], n = 3):
    """"Function to return probable sample name from a path, it recurselvy goes through the path and returns the first element not in the black list."""
    if n == 0:
        return ""

    tmp = path.basename(file_path)
    _d = path.dirname(file_path)

    if all(entry not in tmp for entry in black_list):
        return tmp
    else:
        res = get_sample_name(_d, black_list, n-1)
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

    adata = sc.read_mtx(data_path + 'count_matrix.mtx')

    # reading in gene and cell data
    gene_data = pd.read_csv(data_path + 'all_genes.csv')
    cell_meta = pd.read_csv(data_path + 'cell_metadata.csv')

    # find genes with nan values and filter
    gene_data = gene_data[gene_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()

    # remove genes with nan values and assign gene names
    adata = adata[:,notNa]
    adata.var = gene_data
    adata.var.set_index('gene_name', inplace=True)
    adata.var.index.name = None
    adata.var_names_make_unique()

    # add cell meta data to anndata object
    adata.obs = cell_meta
    adata.obs.set_index('bc_wells', inplace=True)
    adata.obs.index.name = None
    adata.obs_names_make_unique()

    return adata



def human2mouse(genes: List[str]) -> List[str]:

    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
        json={
            'organism':'hsapiens',
            'target':'mmusculus',
            'query':genes,
        }
        )
    df = pd.DataFrame(r.json()['result'], )
    return df.name.replace("N/A", pd.NA).dropna().to_list()





## Technology components

inputs: Dict[str, List|Callable] = {
          "10x":{
                 "files": ['features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz'],
                 "black_list": ["filtered_feature_bc", "raw_feature_bc", "count", "outs"],
                 "raw_name": "raw_feature_bc_matrix",
                 "function": sc.read_10x_mtx
                 },

          "ParseBio":{
                    "files": ["all_genes.csv", "cell_metadata.csv", "count_matrix.mtx"],
                    "black_list": ["DGE_filtered", "DGE_unfiltered"],
                    "function": read_parsebio
                    }
          }


qc_features_fac: Dict[str, List[str]] = {"human": {
                         "mito": ["MT-"],
                         "ribo": ["RBS", "RPL"],
                         "hb": ["^HB[^(P)]"]
                         },
               "mouse": {
                        "mito": ["mt"],
                        "ribo": ["Rps", "Rpl"],
                        "hb": ["^Hb[^(p)]"] # Validate this later
               }
                         }


def reduce_outliers(adata: AnnData, variables: Dict[str, List], subset: bool|None = None) -> pd.Series:
    outlier_dict = {}

    if subset:
        for sample in adata.keys():
            for key in variables[sample].keys():
                if key in adata.obs.columns:
                    if len(variables[key]) == 2:
                        try:
                            adata.obs.loc[adata.obs["sample"] == subset, f"{key}_outlier"] = adata.obs.loc[adata.obs["sample"] == subset, key].lt(variables[key][0]) | \
                                adata.obs.loc[adata.obs["sample"] == subset, key].gt(variables[key][1])
                            outlier_dict[key] = adata.obs[f"{key}_outlier"] 
                        except:
                            adata.obs.loc[adata.obs["sample"] == subset, f"{key}_outlier"] = outlier_dict[key]
                            outlier_dict[key] = adata.obs[f"{key}_outlier"] 
                    else:
                        raise ValueError("Provide a list of length 2 for the lower and upper bound of the QC-variable.")
                else:
                    raise KeyError("the provided QC variable does not exist in the data, check the variable names again.")
            
        return reduce(lambda x, y: x or y, zip(outlier_dict.values()))

    for key in variables.keys():
        if key in adata.obs.columns:
            if len(variables[key]) == 2:
                outlier_dict[key] = adata.obs[key].lt(variables[key][0]) | adata.obs[key].gt(variables[key][1])
                adata.obs[f"{key}_outlier"] = outlier_dict[key]
            else:
                raise ValueError("Provide a list of length 2 for the lower and upper bound of the QC-variable.")
        else:
            raise KeyError("the provided QC variable does not exist in the data, check the variable names again.")

    return reduce(lambda x, y: x or y, zip(outlier_dict.values()))
