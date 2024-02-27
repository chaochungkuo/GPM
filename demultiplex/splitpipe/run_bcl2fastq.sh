# screen -S fastq
# bash run_splitpipe.sh

############################################################
# Demultiplex the reads for sublibraries
PATH_BCL2FASTQ \
  --no-lane-splitting \
  --runfolder-dir PROJECT_BCL_PATH \
  --output-dir ./sublibraries_FASTQ \
  --sample-sheet ./samplesheet_bcl2fastq.csv \
  --processing-threads N_CORES

###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j N_CORES "fastqc {} -o ./fastqc"

###### Running fastq_screen ###############################
mkdir -p ./fastq_screen
fastq_screen --outdir ./fastq_screen --threads N_CORES */*.fastq.gz *.fastq.gz

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f ./fastq_screen/ . -o ./multiqc