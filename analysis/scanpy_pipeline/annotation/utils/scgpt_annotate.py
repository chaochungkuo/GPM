import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1" 
os.environ["WORLD_SIZE"] = "1"
import logging 
logging.basicConfig(level=logging.ERROR)

import torch
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mode
import scanpy as sc
import sklearn
import warnings
import faiss
import scgpt as scg
import argparse
from utils.build_atlas_index_faiss import load_index, vote


###############################################################################################################
################################################ Argument Parsing ################################################
###############################################################################################################


parser = argparse.ArgumentParser(
                    prog='scGPT_annotate',
                    description='Annotates an input scanpy object with scGPT CT model and saves the results to scGPT columns',
                    epilog='')


parser.add_argument('-i', '--input', dest='filename',required=True,  type=str, help='input path of scanpy file on disk.')           
parser.add_argument('-m', '--model-path', dest= 'model_path', required=True, type=str, help='path for scGPT model')
parser.add_argument('-n', '--index', dest= 'index_path', required=True, type=str, help='path for scGPT model')
parser.add_argument('-c', '--column',  dest= 'column_name', default='scGPT', type=str,  help='column to save the results to (default: scGPT)')
parser.add_argument('-o', '--output', type=str, dest='ouput_path', default='')

args = parser.parse_args()

if len(args.ouput_path) == 0:
    output_path = args.filename
else:
    output_path = args.ouput_path


# print("outputpath: ",output_path,"column_name:", args.column_name, "filename: ", args.filename,  "model_path", args.model_path)


model_dir = Path(args.model_path)
adata = sc.read_h5ad(args.filename)
cell_type_key = "scGPT_metainfo"
gene_col = "index"
adata.obs[cell_type_key] = 0



###############################################################################################################
################################################ Gene Embeding ################################################
###############################################################################################################
try:
    adata_embed = scg.tasks.embed_data(
        adata,
        model_dir,
        gene_col=gene_col,
        obs_to_save=cell_type_key,  # optional arg, only for saving metainfo
        batch_size=1024,
        return_new_adata=True,
    )
    adata_embed = adata_embed.X
except:
    exit(code=1)


###############################################################################################################
################################################ Loading Index ################################################
###############################################################################################################

try:
    use_gpu = faiss.get_num_gpus() > 0
    index, meta_labels = load_index(
        index_dir=args.index_path,
        use_config_file=False,
        use_gpu=True,
    )
    print(f"Loaded index with {index.ntotal} cells")
except:
    exit(code=1)


###############################################################################################################
################################################ Index Search #################################################
###############################################################################################################


k = 50
# test with the first 100 cells
distances, idx = index.search(adata_embed, k)

predict_labels = meta_labels[idx]
from scipy.stats import mode
from tqdm import tqdm

voting = []
for preds in tqdm(predict_labels):
    voting.append(vote(preds, return_prob=False)[0])
voting = np.array(voting)
adata.obs[args.column_name] = voting


#################################################################################################################
################################################ Writing results ################################################
#################################################################################################################

adata.write_h5ad(filename=output_path)