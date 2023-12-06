gpm init --help
gpm init --name YYMMDD_Name1_Name2_Institute_Application --fastq FASTQ_Path --processing nfcore_3mRNAseq
cd YYMMDD_Name1_Name2_Institute_Application/nfcore_3mRNAseq/
bash run_nfcore_3mrnaseq.sh --test
# gpm init --name Application --fastq FASTQ_Path --processing nfcore_3mRNAseq