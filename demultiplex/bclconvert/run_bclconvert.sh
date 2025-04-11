# Please make sure bcl-convert is available in your environment
# You should review the paramenters and run:
# screen -S bclconvert
# bash run_bclconvert.sh

# Illumina documentation: https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm
# 10X documentation: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing-bcl-convert

# Define the boolean variable for running fastq_screen
run_fastq_screen=true  # Set to true or false as needed

# Get total logical cores and reserve 10 (minimum 1)
TOTAL_CORES=$(nproc)
AVAILABLE_CORES=$(( TOTAL_CORES - 10 > 0 ? TOTAL_CORES - 10 : 1 ))
# Get CPU idle percentage (average over 1 second)
IDLE_PCT=$(mpstat 1 1 | awk '/Average:/ && $12 ~ /[0-9.]+/ { print $12 }')
# Calculate estimated idle CPUs (rounded down)
IDLE_CPUS=$(awk -v idle_pct="$IDLE_PCT" -v cores="$AVAILABLE_CORES" \
    'BEGIN { print int(idle_pct * cores / 100) }')
echo "Total cores      : $TOTAL_CORES"
echo "Usable cores     : $AVAILABLE_CORES (total - 10)"
echo "Idle %           : $IDLE_PCT"
echo "Estimated idle CPUs : $IDLE_CPUS"
# Calculate 1/3 of NJOBS (rounded down)
THIRDJOBS=$(( IDLE_CPUS / 3 ))

# Please execute this command in the directory OUTPUT_DIR
bcl-convert --bcl-input-directory PROJECT_BCL_PATH \
  --output-directory output \
  --sample-sheet ./samplesheet_bclconvert.csv \
  --bcl-sampleproject-subdirectories true \
  --no-lane-splitting true \
  --bcl-num-conversion-threads $THIRDJOBS \
  --bcl-num-compression-threads $THIRDJOBS \
  --bcl-num-decompression-threads $THIRDJOBS


###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 2 -name "*.fastq.gz" | parallel -j $IDLE_CPUS "fastqc {} -o ./fastqc"

###### Running fastq_screen (Conditional Execution) ########
if [ "$run_fastq_screen" = true ]; then
  echo "Running fastq_screen..."
  mkdir -p ./fastq_screen
  find * -maxdepth 2 -name "*.fastq.gz" | parallel -j $IDLE_CPUS "fastq_screen --outdir ./fastq_screen {}"
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