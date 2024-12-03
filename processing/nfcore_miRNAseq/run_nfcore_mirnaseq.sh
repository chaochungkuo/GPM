################## test script ##################
# PATH_NEXTFLOW run nf-core/rnaseq -profile test,docker

############### nf-core pipeline ###############################

PATH_NEXTFLOW run nf-core/smrnaseq -r 2.3.1 \
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
# --mirtrace_species 'mmu'

# If the run with QIAseq miRNA Library QC Spike- In Sequences, please add the following options:
# --mirna_gtf /data/genomes/spikein/QIASeq_miRNAseq_SpikeIn/hsa.gff3 \
# --mature /data/genomes/spikein/QIASeq_miRNAseq_SpikeIn/mature_with_qiaseq_spikein.fa \
# --hairpin /data/genomes/spikein/QIASeq_miRNAseq_SpikeIn/hairpin_with_qiaseq_spikein.fa \