################## test script #########################
# PATH_NEXTFLOW run nf-core/mag -profile test,docker

PATH_NEXTFLOW run nf-core/mag -r 2.4.0 -profile docker --input samplesheet.csv --outdir results \
--skip_prodigal --skip_prokka --skip_metaeuk --skip_binning --skip_binqc