# Please make sure bcl2fastq V.2 is available in your environment
# You should review the paramenters and run:
# screen -S bcl2fastq
# bash run_bcl2fastq.sh

# Please execute this command in the directory OUTPUT_DIR
PATH_BCL2FASTQ \
  --no-lane-splitting \
  --runfolder-dir BCL_PATH \
  --output-dir . \
  --interop-dir ./InterOp/ \
  --sample-sheet ./samplesheet_bcl2fastq.csv \
  --processing-threads N_CORES

###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j N_CORES "fastqc {} -o ./fastqc"

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f . -o ./multiqc