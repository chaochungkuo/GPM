from os import walk, path, mkdir, listdir
from scipy.stats import median_abs_deviation
import pandas as pd
import numpy as np
from anndata import AnnData
from typing import List, Dict, Callable
import scanpy as sc
import requests
from functools import reduce
from numbers import Number
from scipy.stats import median_abs_deviation


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


###############################################################################################################################################################
#######                                                              Outlier Detection                                                                  #######
###############################################################################################################################################################

def reduce_outliers(adata: AnnData, variables: Dict[str, List|Number]|Dict[str, Dict[str, List|Number]], subset: bool = False) -> pd.Series:
    if subset:
        adata.obs["outlier"] = False
        for sample in variables.keys():
            _compute_outlier_sample(adata, variables, sample)
    else:
        _compute_outlier_all(adata, variables)

def _compute_outlier_all(adata: AnnData, variables: Dict[str, List]) -> pd.Series:
    for key in variables.keys():
        if not key in adata.obs.columns:
            raise KeyError("the provided QC variable does not exist in the data, check the variable names again.")


        if isinstance(variables[key], list):
            if len(variables[key]) != 2:
                raise ValueError("Provide a list of length 2 for the lower and upper bound of the QC-variable.")
            min_val = variables[key][0]
            max_val = variables[key][1]
        
        if isinstance(variables[key], Number):
            if not variables[key] > 0:
                raise ValueError("Please provide a positive number of nmads")
            
            min_val = np.median(adata.obs[key]) - (median_abs_deviation(adata.obs[key]) * variables[key])
            max_val = np.median(adata.obs[key]) + (median_abs_deviation(adata.obs[key]) * variables[key])

        adata.obs[f"{key}_outlier"] = adata.obs[key].lt(min_val) | adata.obs[key].gt(max_val)
        adata.obs[f"{key}_outlier"].fillna(False)
        
    adata.obs["outlier"] =  adata.obs[[f"{x}_outlier" for x in variables.keys()]].any(axis = 1)


def _compute_outlier_sample(adata: AnnData, variables: Dict[str, List], sample):
    sample_dict = variables[sample]

    for key in sample_dict.keys():

        if not key in adata.obs.columns:
            raise KeyError("the provided QC variable does not exist in the data, check the variable names again.")

        if not f"{key}_outlier" in adata.obs.columns:
            adata.obs[f"{key}_outlier"] = False

        if isinstance(sample_dict[key], list):
            if len(sample_dict[key]) != 2:
                raise ValueError("Provide a list of length 2 for the lower and upper bound of the QC-variable.")
            
            min_val = sample_dict[key][0]
            max_val = sample_dict[key][1]

        if isinstance(sample_dict[key], Number):
            if not variables[key] > 0:
                raise ValueError("Please provide a positive number of nmads")
            min_val = np.median(adata.obs[key]) - (median_abs_deviation(adata.obs[key]) * variables[key])
            max_val = np.median(adata.obs[key]) + (median_abs_deviation(adata.obs[key]) * variables[key])

        sample_slice =  adata.obs.loc[adata.obs["sample"] == sample, key].lt(min_val) | \
                        adata.obs.loc[adata.obs["sample"] == sample, key].gt(max_val)
        
        adata.obs.loc[adata.obs["sample"] == sample, f"{key}_outlier"] = sample_slice
        adata.obs.loc[adata.obs["sample"] == sample, f"{key}_outlier"].fillna(False)
    
    adata.obs.loc[adata.obs["sample"] == sample, "outlier"] = \
        adata.obs.loc[adata.obs["sample"] == sample, [f"{x}_outlier" for x in sample_dict.keys()]].any(axis = 1)



def get_keys(qc_dict):
    keys_list = []
    if len(qc_dict) > 0 and all(map(lambda x: isinstance(x, (list, Number)), qc_dict.values())):
        return list(qc_dict.keys())
    
    if len(qc_dict) > 0 and all(map(lambda x: isinstance(x, (dict, Number)), qc_dict.values())):
        for key in qc_dict.keys():
            keys_list = keys_list + list(qc_dict[key].keys())
        return list(set(keys_list))
    

    return []

##############################################
############ Normalization Functions #########
##############################################

def get_var_features_num(adata: AnnData, variable_features: int | float) -> int:
    detected_gene_nu = len(adata.var_names)
    if variable_features <= 1:
        return int(detected_gene_nu * variable_features)
    else:
        return min(detected_gene_nu, variable_features)

def is_raw_counts(matrix) -> bool:
    from scipy.sparse import issparse

    if issparse(matrix):
        return matrix.count_nonzero() == matrix.astype("uint32").count_nonzero()
    else:
        return np.count_nonzero(matrix) == np.count_nonzero(matrix.astype("uint32"))