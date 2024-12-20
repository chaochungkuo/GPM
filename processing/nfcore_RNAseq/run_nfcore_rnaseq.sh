#!/bin/bash

nfcore_version=3.17.0

# Check if the --test option is provided
if [ "$#" -gt 0 ] && [ "$1" = "--test" ]; then
  ################## test script #########################
  PATH_NEXTFLOW run nf-core/rnaseq -r $nfcore_version \
  -profile test,docker --outdir results
else
     ################## GPM samplesheet #####################
     # gpm samplesheet-rnaseq -st 'reverse' -sn true -si 2 samplesheet.csv PROJECT_FASTQ_PATH

     ################## Run nfcore pipeline #################

     PATH_NEXTFLOW run nf-core/rnaseq -r ${nfcore_version} -profile docker -c nextflow.config \
          --input samplesheet.csv --outdir results \
          --genome GENCODE_GRCh38_v46  \
          --gencode --featurecounts_group_type gene_type
          # --with_umi --umitools_extract_method "regex" --umitools_bc_pattern2 "^(?P<umi_1>.{8})(?P<discard_1>.{6}).*" # Takara Bio SMARTer® Stranded Total RNA-Seq Kit v3	
fi

###### For rerun the pipeline #################################
# use: -resume

###### Options for genome: ####################################
# GENCODE_GRCh38_v46, GENCODE_GRCh38_v46_ERCC, 
# GENCODE_GRCm39_v35, GENCODE_GRCm39_v35_ERCC
# --gencode --featurecounts_group_type gene_type

###### For runs with ERCC spike-ins ############################
# --genome GENCODE_GRCh38_v46_ERCC --skip_biotype_qc