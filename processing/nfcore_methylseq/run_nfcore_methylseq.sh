#!/bin/bash

nfcore_version=2.6.0

# Check if the --test option is provided
if [ "$#" -gt 0 ] && [ "$1" = "--test" ]; then
  ################## test script #########################
  PATH_NEXTFLOW run nf-core/methylseq -r $nfcore_version \
  -profile test,docker --outdir results
else
  ################## GPM samplesheet #####################
  # Because the format is the same as nfcore scRNAseq, we can use the same
  # command gpm samplesheet-scrnaseq to generate samplesheet.csv
  # gpm samplesheet-scrnaseq -sn true -si 1 samplesheet.csv PROJECT_FASTQ_PATH

  ################## Run nfcore pipeline #################
  PATH_NEXTFLOW run nf-core/methylseq -r $nfcore_version -profile docker -c nextflow.config \
  --input samplesheet.csv --outdir results --genome GENCODE_GRCh38_v44

fi

###### For rerun the pipeline #################################
# use: -resume

###### Options for genome: ####################################
# GENCODE_GRCh38_v44, GENCODE_GRCh38_v44_ERCC, 
# GENCODE_GRCm39_v33, GENCODE_GRCm39_v33_ERCC
# --gencode --featurecounts_group_type gene_type
