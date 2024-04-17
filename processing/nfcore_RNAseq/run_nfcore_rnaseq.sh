#!/bin/bash

nfcore_version=3.12.0

# Check if the --test option is provided
if [ "$#" -gt 0 ] && [ "$1" = "--test" ]; then
  ################## test script #########################
  PATH_NEXTFLOW run nf-core/rnaseq -r $nfcore_version \
  -profile test,docker --outdir results
else
     ################## GPM samplesheet #####################
     # gpm samplesheet-rnaseq -st 'reverse' -sn true -si 2 samplesheet.csv PROJECT_FASTQ_PATH

     ################## Run nfcore pipeline #################

     PATH_NEXTFLOW run nf-core/rnaseq -r 3.12.0 -profile docker -c nextflow.config \
          --input samplesheet.csv --outdir results \
          --genome GENCODE_GRCh38_v44  \
          --gencode --featurecounts_group_type gene_type
fi

###### For rerun the pipeline #################################
# use: -resume

###### Options for genome: ####################################
# GENCODE_GRCh38_v44, GENCODE_GRCh38_v44_ERCC, 
# GENCODE_GRCm39_v33, GENCODE_GRCm39_v33_ERCC
# --gencode --featurecounts_group_type gene_type

###### For runs with ERCC spike-ins ############################
# --genome GENCODE_GRCh38_v44_ERCC --skip_biotype_qc