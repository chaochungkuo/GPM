# Usage example

```bash
Rscript ./sc_conversion.R -i /path/to/input.{h5ad, RDS} -o /path/to/output.{h5ad, RDS} -t {seurat, sce, loupe} 
```

Not all the conversion pathways is working

# Implementation table

| Input format               | Output format              | Supported                                          |
| -------------------------- | -------------------------- | -------------------------------------------------- |
| AnnData (.h5ad)            | SingleCellExperiment (RDS) | Yes :white_check_mark:                             |
| SingleCellExperiment (RDS) | AnnData (.h5ad)            | Yes :white_check_mark:                             |
| Seurat (RDS)               | AnnData (.h5ad)            | Partial (Bugs can occur ) :heavy_exclamation_mark: |
| AnnData (.h5ad)            | Seurat (RDS)               | NO :x:                                             |
| Seurat (RDS)               | SingleCellExperiment (RDS) | NO :x:                                             |
| Any                        | Loupe                      | NO :x:                                             |
