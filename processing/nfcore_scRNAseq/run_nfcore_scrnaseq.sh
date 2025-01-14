################## test script ##################
# PATH_NEXTFLOW run nf-core/scrnaseq -profile test,docker

################## GPM samplesheet #####################
# gpm samplesheet-scrnaseq -sn true -si 2 samplesheet.csv PROJECT_FASTQ_PATH

PATH_NEXTFLOW run nf-core/scrnaseq -r 2.5.1 -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --cellranger_index REFDATA_CELLRANGER/refdata-gex-GRCh38-2020-A \
    --aligner cellranger  \
    --protocol 10XV3 \
    --skip_emptydrops \ # Skip CellBender due to long running time on CPU
    --multiqc_title PROJECT_PROJECT_NAME


# Options for --genome:
# GRCh38, GRCm39, mm10
# Example for --cellranger_index: (Specify a pre-calculated cellranger index. Has to correspond and match the genome parameter's type)
# '/data/shared_env/10xGenomics/refdata-gex-GRCh38-2020-A'
# '/data/shared_env/10xGenomics/refdata-gex-mm10-2020-A'
# '/data/shared_env/10xGenomics/refdata-gex-GRCh38-and-mm10-2020-A'

