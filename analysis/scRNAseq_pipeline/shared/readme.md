
# sc_convert
sc_convert is a script for disk-based conversion of scRNA-seq data between Seurat, SingleCellExperiment, and AnnData formats. It also supports conversion to Loupe format.
All conversions from or to Seurat, SingleCellExperiment, and AnnData are supported. Loupe is supported only as a final output.

## Installation

`sc_convert` uses `docker` to conterarize all tools required for running the conversion. To install `sc_convert`, you need to have `docker` and `docker compose` installed on your machine.

```sh
git clone https://github.com/chaochungkuo/GPM
cd GPM/analysis/scanpy_pipeline/shared
chmod +x sc_convert.sh
docker compose build
```


## Usage example

```sh
bash sc_convert.sh --input|-i /path/to/input.{h5ad, RDS}  --ouput|-o /path/to/output.{h5ad, RDS, cloupe} --from|-f {seurat, sce, anndata}  -t {seurat, sce, loupe, anndata} 
```





