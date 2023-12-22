gpm init --help
gpm init --name 231221_Name1_Name2_Institute_3mRNAseq --fastq FASTQ_Path --processing nfcore_3mRNAseq

cd 231221_Name1_Name2_Institute_3mRNAseq
gpm processing -fq FASTQ_PATH -p nfcore_RNAseq project.ini
gpm processing -fq FASTQ_PATH -p nfcore_miRNAseq project.ini

# Analysis
gpm analysis --list project.ini
gpm analysis --report RNAseq project.ini
gpm analysis --add DGEA_RNAseq project.ini
gpm analysis --add GO_analysis project.ini
gpm analysis --add GSEA project.ini
cd ..
# bash run_nfcore_3mrnaseq.sh --test
# gpm init --name Application --fastq FASTQ_Path --processing nfcore_3mRNAseq