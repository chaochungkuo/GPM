
# Quick start 

You should have the following installed:
- [`GPM`](https://github.com/chaochungkuo/GPM) 
- [`Pixi`](https://pixi.sh/latest/)
  
## Using GPM as project manager.

**Run all the preprocessing steps at once:**

```sh
gpm analysis <path/to/project.ini> --add preprocessing --report scRNAseq  
cd <path/to/project>/analysis/scRNAseq_pipeline
pixi install
pixi run render_all
```

**to render the notebooks one by one:**
```sh
pixi run render <notebook_name>
```

# scRNAseq_pipeline
a semi-supervised pipeline for single cell RNA sequencing data analysis. It is designed to be modular and to automate most of scRNA-seq preprocessing tasks. It has a collection of Jupyter notebooks that that mix the code and visualization of the results. Most of notebook comes with sensible parameters that should generate reasonable results for most of the datasets. For manual tuning, the user can modify the parameters in the notebook and re-run the analysis.

## Components
- `Config.toml `: holds project-wide information as organism, technology and the location of the data.
- `modules`: each module have one or more jupyter notebook with the same structure. 
  - The first code cell contain the parameter for running the analysis with explanation of each parameter. The rest of the notebook is the code that runs the analysis.
  - The parameters can be freely adjusted and the notebook can be ran all at once or cell by cell.
  - Each module is isolated as a pixi environment and can be run independently.
- `resources`: contains the data that is used in the analysis as cell-cycle genes and stores ML models if used.
- `save`: The primary location for storing AnnData objects that are re-used through out the analysis and other results as Differential expression tables.

## Modules
### Preprocessing
- **QC**: Quality control module that filters out low quality cells and genes. It also detects doublets and ambient RNA.
- **Normalization & Clustering**: Normalizes the data by library size and log transform the data.
- **Annotation**: Annotates the cells with cell type labels.
### Downstream analysis 
Can be only ran after preprocessing
- **Differential_expression**: Performs differential expression analysis between two groups of cells.
- **RNA_velocity**: Estimates the RNA velocity of the cells.  


# TODO
- [ ]  Add CI/CD tests based on scanpy databases.
  

TODO - QC:
- [X] Make QC per sample by default.
- [X] Simpify the plotting of the QC module and clean up the manifest file.
- [X] Turn on doublet detection and Ambient RNA detection.
- [X] Add Documentation for the usage of the QC module.
- [ ] Refactor the sample discovery code


 