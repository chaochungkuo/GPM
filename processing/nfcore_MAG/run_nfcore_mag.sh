################## test script #########################
# nextflow run nf-core/mag -profile test,docker

nextflow run nf-core/mag -r 2.4.0 -profile docker --input samplesheet.csv --outdir results \
--skip_prodigal --skip_prokka --skip_metaeuk --skip_binning --skip_binqc