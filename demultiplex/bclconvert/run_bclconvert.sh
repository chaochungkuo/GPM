# Please make sure bcl-convert is available in your environment
# You should review the paramenters and run:
# screen -S bclconvert
# bash run_bclconvert.sh

# Define the boolean variable for running fastq_screen
run_fastq_screen=true  # Set to true or false as needed

# Get the nummber of jobs as ncors - 10 or 1
NJOBS=$(( $(nproc) - 10 > 0 ? $(nproc) - 10 : 1 ))

# Please execute this command in the directory OUTPUT_DIR
bcl-convert -bcl-input-directory PROJECT_BCL_PATH \
  --output-directory . \
  --sample-sheet ./samplesheet_convert.csv \
  --bcl-sampleproject-subdirectories true \
  --sample-name-column-enabled true \
  --no-lane-splitting true \
  --bcl-num-conversion-threads $NJOBS \
  --bcl-num-compression-threads $NJOBS \
  --bcl-num-decompression-threads $NJOBS


###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j $NJOBS "fastqc {} -o ./fastqc"

###### Running fastq_screen (Conditional Execution) ########
if [ "$run_fastq_screen" = true ]; then
  echo "Running fastq_screen..."
  mkdir -p ./fastq_screen
  find * -maxdepth 1 -name "*.fastq.gz" | parallel -j $NJOBS "fastq_screen --outdir ./fastq_screen {}"
else
  echo "Skipping fastq_screen as run_fastq_screen is set to false."
fi

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f . -o ./multiqc

# Cleanup
rm -r ./fastqc
if [ -d "./fastq_screen" ]; then
  rm -r ./fastq_screen
fi