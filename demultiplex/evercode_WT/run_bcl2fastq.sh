# screen -S fastq
# bash run_splitpipe.sh

############################################################

# Get the nummber of jobs as ncors - 10 or 1
NJOBS=$(( $(nproc) - 10 > 0 ? $(nproc) - 10 : 1 ))


# Demultiplex the reads for sublibraries
PATH_BCL2FASTQ \
  --no-lane-splitting \
  --runfolder-dir PROJECT_BCL_PATH \
  --output-dir ./sublibraries_FASTQ \
  --sample-sheet ./samplesheet_bcl2fastq.csv \
  --processing-threads $NJOBS

###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j $NJOBS "fastqc {} -o ./fastqc"

###### Running fastq_screen ###############################
mkdir -p ./fastq_screen
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j $NJOBS "fastq_screen --outdir ./fastq_screen {}"

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f ./fastq_screen/ . -o ./multiqc

rm -r ./fastqc ./fastq_screen