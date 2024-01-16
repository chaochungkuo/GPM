conda activate schashtag

EXP_CELL=5000
NUM_THREADS=N_CORES
FASTQ_DIR="PROJECT_FASTQ_PATH"
for FASTQ in ${FASTQ_DIR}/*R1_001.fastq.gz; do
  name=$(basename "$FASTQ" | cut -d'_' -f1,2)
  echo $name
  R1=$(find "$FASTQ_DIR" -type f -name ${name}*R1*.fastq.gz)
  R2=$(find "$FASTQ_DIR" -type f -name ${name}*R2*.fastq.gz)
  # echo $R1
  # echo $R2
  CITE-seq-Count -R1 $R1 -R2 $R2 -t tag.csv -T $NUM_THREADS \
  -cbf 1 -cbl 16 -umif 17 -umil 28 -cells $EXP_CELL -o $name
done

# For 10X Genomics: -cbf 1 -cbl 16 -umif 17 -umil 28
# https://citeseq.files.wordpress.com/2019/02/cite-seq_and_hashing_protocol_190213.pdf