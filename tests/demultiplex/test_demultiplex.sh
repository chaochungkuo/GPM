rm ~/gpmdata/*
pip install ../../
# rm -r /Users/jovesus/Downloads/GPM
mkdir YYMMDD_M05588_0232_000000000-CLL8M
tar -xzf YYMMDD_M05588_0232_000000000-CLL8M.tgz -C YYMMDD_M05588_0232_000000000-CLL8M
gpm demultiplex --raw YYMMDD_M05588_0232_000000000-CLL8M --output ~/demultiplex_output/ --method bcl2fastq
