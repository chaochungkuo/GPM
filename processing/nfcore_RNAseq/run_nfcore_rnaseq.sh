################## test script ##################
# nextflow run nf-core/rnaseq -r 3.12.0 -profile test,docker --outdir results

################## GPM samplesheet #####################
# gpm samplesheet -st 'auto' -sn true -si  samplesheet.csv FASTQ_DIR

nextflow run nf-core/rnaseq -r 3.12.0 -profile docker -c nextflow.config \
     --input samplesheet.csv --outdir results \
     --genome GENCODE_GRCh38_v44  \
     --max_cpus 45 --gencode --featurecounts_group_type gene_type

# For rerun the pipeline, use: -resume
# Options for genome: GENCODE_GRCh38_v44, GENCODE_GRCh38_v44_ERCC