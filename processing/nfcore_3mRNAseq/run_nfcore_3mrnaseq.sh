################## test script #########################
# nextflow run nf-core/rnaseq -r 3.12.0 -profile test,docker --outdir results

################## GPM samplesheet #####################
# gpm samplesheet -st 'forward' -sn true -si samplesheet.csv FASTQ_DIR

nextflow run nf-core/rnaseq -r 3.12.0 -profile docker -c nextflow.config \
     --input samplesheet.csv --outdir results \
     --genome GENCODE_GRCh38_v44 --gencode --featurecounts_group_type gene_type \
     --max_cpus 45 \
     --extra_salmon_quant_args="--noLengthCorrection" \
     --extra_star_align_args="--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA" \
     --with_umi --umitools_extract_method="regex" --umitools_bc_pattern="^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"

# For rerun the pipeline, use: -resume

# Options for genome: GENCODE_GRCh38_v44, GENCODE_GRCh38_v44_ERCC

# QuantSeq 3’mRNA-Seq V2 Library Prep Kit with UDI:
# --extra_star_align_args="--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA" \
# --with_umi --umitools_extract_method="regex" --umitools_bc_pattern="^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"

# For any 3'mRNAseq: --extra_salmon_quant_args="--noLengthCorrection"

# For runs with ERCC spike-ins:
# --genome GENCODE_GRCh38_v44_ERCC --skip_biotype_qc

# --gencode --featurecounts_group_type gene_type