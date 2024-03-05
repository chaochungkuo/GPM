################## test script ##################
# PATH_NEXTFLOW run nf-core/rnaseq -profile test,docker


######################################################## Indpendent UMI extraction ########################################################
# for FASTQ in PROJECT_FASTQ_PATH/*.fastq.gz
# mkdir -p PROJECT_FASTQ_PATH/UMI_trimmed
# ls PROJECT_FASTQ_PATH/*.fastq.gz | parallel -j 40 'filename=$(basename -- {});\
# echo $filename; TRIMMED_FASTQ="PROJECT_FASTQ_PATH/UMI_trimmed/${filename}"; \
# echo $TRIMMED_FASTQ; umi_tools extract --stdin={} --stdout=${TRIMMED_FASTQ} \
# --extract-method=regex --bc-pattern=".+AACTGTAGGCACCATCAAT{s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)"'



################## GPM samplesheet #####################
# gpm samplesheet -st 'forward' samplesheet.csv PROJECT_FASTQ_PATH/UMI_trimmed

######################## nf-core pipeline #####################################
# Run the following in the PROJECT_FASTQ_PATH
# nextflow run nf-core/smrnaseq -r 2.3.0 \
#      -profile docker \
#      --input samplesheet.csv \
#      --genome 'GRCh38' \
#      --mirtrace_species 'hsa' \
#      --protocol 'qiaseq' \
#      --outdir results \
#      --save_reference \


######################################################## Combined Pipeline #################################################################

################## GPM samplesheet ############################
# gpm samplesheet -st 'forward' samplesheet.csv PROJECT_FASTQ_PATH


############### nf-core pipeline ############################### 
PATH_NEXTFLOW run nf-core/smrnaseq -r 2.3.0 \
     -profile docker \
     --input samplesheet.csv \
     --genome 'GRCh38' \
     --mirtrace_species 'hsa' \
     --protocol 'qiaseq' \
     --outdir results \
     --save_reference \
     --with_umi \
     --umitools_extract_method regex \
     --umitools_bc_pattern '.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)' # Regex pattern for Qiagen_QIAseq_miRNA

# Options for --genome:
# gencode_hg38, gencode_mm10, hg38, mm10
# --mirna_gtf /data/genomes/hg38/miRNA/hsa.gff3 --mature /data/genomes/hg38/miRNA/mature.fa.gz --hairpin /data/genomes/hg38/miRNA/hairpin.fa.gz