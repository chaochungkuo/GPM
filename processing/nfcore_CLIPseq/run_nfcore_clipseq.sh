################## test script ##################
# nextflow run nf-core/clipseq -profile test, docker

nextflow run goodwright/clipseq -latest -profile docker \
     --samplesheet samples.csv --outdir results \
     --genome GENCODE_GRCm38_v25 --fasta '/data/genomes/GRCm38/GRCm38.p6.genome.fa' \
     --gtf '/data/genomes/GRCm38/gencode.vM25.annotation.gtf' \
     --smrna_fasta Mus_musculus.smallRNA.fa.gz \
     --peakcaller "icount,paraclu,pureclip,piranha" --skip_umi_dedupe true --motif

# For rerun the pipeline, use: -resume
# --smrna_fasta, search FASTQ files here:
# https://github.com/nf-core/clipseq/tree/master/assets/small_rna