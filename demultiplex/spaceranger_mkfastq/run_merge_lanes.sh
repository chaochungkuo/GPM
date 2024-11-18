#!/bin/bash
# This script will:
#  1) Merge the lanes by samples
#  2) Run FASTQC on the merged FASTQs
#  3) Run MultiQC
#  4) Delete merged FASTQs (Because these files are only used for QC, not for downstream analyses)
OUT_DIR="$1"
MQC_DIR="$2"
FASTQ_DIR=$(find ${OUT_DIR}/outs/fastq_path/ -maxdepth 1 -type d -name "H*")
samples=$(find ${OUT_DIR}/outs/fastq_path/ -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/_S[0-9]_.*//' | sort | uniq)
MERGED_DIR="./merged_fastq"
###### Merging lanes #######################################
mkdir -p $MERGED_DIR
echo $samples
for sample in $samples; do
  if [[ "$sample" != "Undetermined" ]]; then
    echo -e "${sample}\tMerging R1"
    cat ${FASTQ_DIR}/${sample}_S*_L*_R1_001.fastq.gz > ${MERGED_DIR}/${sample}_Merged_R1_001.fastq.gz
    echo -e "${sample}\tMerging R2"
    cat ${FASTQ_DIR}/${sample}_S*_L*_R2_001.fastq.gz > ${MERGED_DIR}/${sample}_Merged_R2_001.fastq.gz
  else
    echo -e "Undetermined\tMerging R1"
    cat ${OUT_DIR}/outs/fastq_path/Undetermined_S*_L*_R1_001.fastq.gz > ${MERGED_DIR}/Undetermined_Merged_R1_001.fastq.gz
    echo -e "Undetermined\tMerging R2"
    cat ${OUT_DIR}/outs/fastq_path/Undetermined_S*_L*_R2_001.fastq.gz > ${MERGED_DIR}/Undetermined_Merged_R2_001.fastq.gz
  fi
done

###### Running FASTQC ######################################
mkdir -p ./fastqc
find $MERGED_DIR -maxdepth 1 -name "*.fastq.gz" | parallel -j 10 "fastqc {} -o ./fastqc"

###### Running fastq_screen ######################################
mkdir -p ./fastq_screen
fastq_screen --outdir ./fastq_screen --threads 10 ${MERGED_DIR}/*.fastq.gz

###### Running MultiQC #####################################
mkdir -p ${MQC_DIR}
multiqc -f . ./fastq_screen/ -o ./${MQC_DIR}

rm -r $MERGED_DIR
rm -r ./fastqc ./fastq_screen