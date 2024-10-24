# Version 0.2.1
## QC
- Added `ScaleBio` to `AutoDiscover` module to automatically discover samples from simple `ScaleBio` nextflow pipeline output.
# Version 0.2.0
## QC
- Changed the `CONTACT_SAMPLES` variable to `QC_PER_SAMPLE`. \n  now all samples are concatenated and the `QC_PER_SAMPLE` controls if the outlier cells are detected in aggregate or per sample. 
- Two types of QC variables where added, `MAX_ONLY`and `LOG_TRANSFORMED`. `MAX_ONLY` Variables as mitchondrial and ribosomal genes where only the upper threshold is applied when using number of MADS as threshold.  `LOG_TRANSFORMED` where the log transformed values are used for the outlier detection as `total_counts` and `n_genes_by_counts`.
- Variable regression was removed from QC module. It is not considered part of the standard QC process to regress variables.
- reading of `10X` output was chnaged to read `.h5` files by default instead of reading separate files for `barcodes`, `features`, and `matrix` .
- Sample Autodiscovery was complelty re-factored and it is not more robust and faster in sample path finding and sample names, it can also detect the output of several output formats from 'CellRanger' as 'Multi' and 'Aggr' outputs. in addition to Singerlon and ParseBio outputs.
- New section `basic.raw_samples` was added to the `config.toml` to specify the raw samples that can be used for Ambient RNA detection and doublet detection.
  
# Version 0.1.0
- Initial version of the workflow.
