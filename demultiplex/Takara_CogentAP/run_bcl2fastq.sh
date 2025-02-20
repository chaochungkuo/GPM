# Please make sure bcl2fastq V.2 is available in your environment
# You should review the paramenters and run:
# screen -S bcl2fastq
# bash run_bcl2fastq.sh

# Define the boolean variable for running fastq_screen
run_fastq_screen=true  # Set to true or false as needed

# Get the nummber of jobs as ncors - 10 or 1
NJOBS=$(( $(nproc) - 10 > 0 ? $(nproc) - 10 : 1 ))

# Please execute this command in the directory OUTPUT_DIR
PATH_BCL2FASTQ \
  --no-lane-splitting \
  --runfolder-dir PROJECT_BCL_PATH \
  --output-dir . \
  --interop-dir ./InterOp/ \
  --sample-sheet ./samplesheet_dummy.csv \
  --processing-threads $NJOBS

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