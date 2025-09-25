# Please make sure bcl-convert is available in your environment
# You should review the paramenters and run:
# screen -S bclconvert
# bash run_bclconvert.sh

# Illumina documentation: https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm
# 10X documentation: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing-bcl-convert

# Define the boolean variable for running fastq_screen
run_fastq_screen=false  # Set to true or false as needed

# Get the nummber of jobs as ncors - 10 or 1
NJOBS=$(( $(nproc) - 10 > 0 ? $(nproc) - 10 : 1 ))

# Please execute this command in the directory OUTPUT_DIR
bcl-convert --bcl-input-directory PROJECT_BCL_PATH \
  --output-directory output \
  --sample-sheet ./samplesheet_bclconvert.csv \
  --bcl-sampleproject-subdirectories false \
  --no-lane-splitting false \
  --bcl-num-conversion-threads 20 \
  --bcl-num-compression-threads 10 \
  --bcl-num-decompression-threads 2


###### Running FASTQC ######################################
mkdir -p ./fastqc
find ./output -maxdepth 1 -name "*.fastq.gz" | parallel -j $NJOBS "fastqc {} -o ./fastqc"

###### Running fastq_screen (Conditional Execution) ########
#if [ "$run_fastq_screen" = true ]; then
#  echo "Running fastq_screen..."
#  mkdir -p ./fastq_screen
#  find * -maxdepth 2 -name "*.fastq.gz" | parallel -j $NJOBS "fastq_screen --outdir ./fastq_screen {}"
#else
#  echo "Skipping fastq_screen as run_fastq_screen is set to false."
#fi

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f . -o ./multiqc

# Cleanup
#rm -r ./fastqc
#if [ -d "./fastq_screen" ]; then
#  rm -r ./fastq_screen
#fi
