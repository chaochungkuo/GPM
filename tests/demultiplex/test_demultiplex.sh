mkdir -p YYMMDD_M05588_0232_000000000-CLL8M
tar -xzf YYMMDD_M05588_0232_000000000-CLL8M.tgz -C .


gpm demultiplex --help
mkdir -p FASTQs
gpm demultiplex --raw YYMMDD_M05588_0232_000000000-CLL8M --output FASTQs --method bcl2fastq
cp -f samples.csv FASTQs/YYMMDD_M05588_0232_000000000-CLL8M/samplesheet_bcl2fastq.csv

cd FASTQs/YYMMDD_M05588_0232_000000000-CLL8M/
bash run_bcl2fastq.sh
cd ../../

rm -r FASTQs
rm -r YYMMDD_M05588_0232_000000000-CLL8M
