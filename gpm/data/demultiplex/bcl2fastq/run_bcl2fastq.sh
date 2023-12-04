# Please make sure bcl2fastq V.2 is available in your environment
# You should review the paramenters and run:
# screen -S bcl2fastq
# bash run_bcl2fastq.sh

# Please execute this command in the directory OUTPUT_DIR
bcl2fastq \
  --no-lane-splitting \
  --runfolder-dir BCL_PATH \
  --output-dir . \
  --interop-dir ./InterOp/ \
  --sample-sheet ./samplesheet_bcl2fastq.csv \
  --processing-threads 30

###### Running FASTQC ######################################
mkdir -p ./fastqc
find * -maxdepth 1 -name "*.fastq.gz" | parallel -j 30 "fastqc {} -o ./fastqc"

###### Running MultiQC #####################################
mkdir -p multiqc
multiqc -f . -o ./multiqc