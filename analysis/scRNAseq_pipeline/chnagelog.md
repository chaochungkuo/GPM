# Version 0.2.0
## QC
- Changed the `CONTACT_SAMPLES` variable to `QC_PER_SAMPLE`.\n now all samples are concatenated and the `QC_PER_SAMPLE` controls if the outlier cells are detected in aggregate or per sample. 
- Two types of QC variables where added, `MAX_ONLY`and `LOG_TRANSFORMED`. `MAX_ONLY` Variables as mitchondrial and ribosomal genes where only the upper threshold is applied when using number of MADS as threshold.  `LOG_TRANSFORMED` where the log transformed values are used for the outlier detection as `total_counts` and `n_genes_by_counts`.
- Variable regression was removed from QC module. It is not considered part of the standard QC process to regress variables.

# Version 0.1.0
- Initial version of the workflow.
