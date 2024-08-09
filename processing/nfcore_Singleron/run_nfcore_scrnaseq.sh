################## test script ##################
# nextflow run singleron-RD/scrna -profile test,docker --outdir results

################## GPM samplesheet #####################
# gpm samplesheet-scrnaseq -sn true -si 2 samplesheet.csv PROJECT_FASTQ_PATH


nextflow run singleron-RD/scrna -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --fasta "/data/genomes/GRCh38/GRCh38.primary_assembly.genome.fa" \
    --gtf "/data/genomes/GRCh38/gencode.v44.basic.annotation.gtf" \
    --genome_name "GENCODE_GRCh38_v44" \
    --run_subsample TRUE \
    --run_fastqc TRUE \
    --publish_dir_mode copy 

    
    

# Options:
# --star_genome "/data/genomes/singleron_genome/GENCODE_GRCh38_v44" 
#--fasta "/data/genomes/GRCh38/GRCh38.primary_assembly.genome.fa" \
# --gtf "/data/genomes/GRCh38/gencode.v44.basic.annotation.gtf" \
#--genome_name "GENCODE_GRCh38_v44" \
