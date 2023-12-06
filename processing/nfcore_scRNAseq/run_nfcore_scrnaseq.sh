################## test script ##################
# nextflow run nf-core/scrnaseq -profile test,docker


nextflow run nf-core/scrnaseq -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --cellranger_index {cellranger-genome-reference} \ 
    --aligner cellranger  \
    --protocol 10XV3


# Options for --genome:
# GRCh38, GRCm39, mm10
# Example for --cellranger_index: (Specify a pre-calculated cellranger index. Has to correspond and match the genome parameter's type)
# '/data/shared_env/10xGenomics/refdata-gex-GRCh38-2020-A'
