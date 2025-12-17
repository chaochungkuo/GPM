uv pip install -e ..
gpm --help

# cd demultiplex
# bash test_demultiplex.sh
# cd ..

# cd processing
# bash test_processing.sh
# cd ..

cd data

gpm init -fq data/sample_fastq -n 241212_Name1_Name2_Test_3mRNAseq -p nfcore_3mRNAseq

cd 241212_Name1_Name2_Test_3mRNAseq/nfcore_3mRNAseq

cd ../../../
