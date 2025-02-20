# This workflow shows the demultiplexing process for Shasta total RNAseq data

##Step1##
# apply run_bcl2fastq.sh to generated 2 large fastq files 


##Step2##
#Demultiplex dry run with cogentAP

cogent rna demux --dry_run -f PROJECT_FASTQ_PATH/Undetermined_S0_R1_001.fastq.gz \
    -p PROJECT_FASTQ_PATH/Undetermined_S0_R2_001.fastq.gz \
    -b /data/shared/takara/CogentAP/config/well_list_shasta_total_rna.csv \
    -t shasta_total_rna \
    -o PROJECT_FASTQ_PATH/dry_run_output

##Step3##
#knee plot to evaluate cut off for the singlecell demultiplex
# apply KneePlot.R
# the calculation of the inflection point is not very good, so it's better to use the interactive plot to make a better decision

##Step4## 
# demultiplex with the cut off 

cogent rna demux --dry_run -f PROJECT_FASTQ_PATH/Undetermined_S0_R1_001.fastq.gz \
    -p PROJECT_FASTQ_PATH/Undetermined_S0_R2_001.fastq.gz \
    -b /data/shared/takara/CogentAP/config/well_list_shasta_total_rna.csv \
    -t shasta_total_rna \
    --use_barcodes 40000 \ #use your own dry run result here
    --min_reads 10000    \ #use your own dry run result here
    -o PROJECT_FASTQ_PATH/output

##Step4##
# map the samples using the sample barcodes
map_samples.py -c COUNTS_FILE -o BARCODE_SAMPLE_MAP_FILE -s SAMPLE_LAYOUT_FILE 
