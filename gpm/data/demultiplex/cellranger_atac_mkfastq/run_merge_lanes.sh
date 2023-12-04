#!/bin/bash
# This script will:
#  1) Merge the lanes by samples
#  2) Run FASTQC on the merged FASTQs
#  3) Run MultiQC
#  4) Delete merged FASTQs (Because these files are only used for QC, not for downstream analyses)

FASTQ_Inout="../CELLRANGER_FASTQ_PATH"
FASTQ_Output="../OUTPUT_DIR"

###### Merging lanes #######################################
mkdir -p $FASTQ_Output

samples=$(ls ${FASTQ_Inout}/*.fastq.gz | xargs basename -a | sed 's/_S[0-9]_.*//' | uniq)
echo $samples

for sample in $samples; do

  echo -e "${sample}\tMerging R1"
  cat ${FASTQ_Inout}/${sample}_S*_L*_R1_001.fastq.gz > ${FASTQ_Output}/${sample}_Merged_R1_001.fastq.gz
  echo -e "${sample}\tMerging R2"
  cat ${FASTQ_Inout}/${sample}_S*_L*_R2_001.fastq.gz > ${FASTQ_Output}/${sample}_Merged_R2_001.fastq.gz

  
done

###### Running FASTQC ######################################
mkdir -p ./fastqc
find $FASTQ_Output -maxdepth 1 -name "*.fastq.gz" | parallel -j 30 "fastqc {} -o ./fastqc"

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f . -o ./multiqc

rm -r $FASTQ_Output