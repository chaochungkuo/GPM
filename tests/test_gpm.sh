pip install ..
gpm --help

cd demultiplex
bash test_demultiplex.sh
cd ..

cd processing
bash test_processing.sh
cd ..

# rm -r 231208_name1_name2_insti_RNAseq/
gpm init -n 231208_name1_name2_insti_RNAseq
cd 231208_name1_name2_insti_RNAseq
gpm processing -fq FASTQ_PATH -p nfcore_RNAseq project.ini
gpm processing -fq FASTQ_PATH -p nfcore_miRNAseq project.ini
# cat project.ini
cd ..