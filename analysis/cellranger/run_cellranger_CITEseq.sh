#!/bin/bash -l

#make sure cellranger-7.1.0 is available in your environment
#otherwise scource shared enviroment from /data/shared_env/shared_paths.sh
#this script can be run directly from shell

fastqDir="/data/fastq/220311_NB501289_0601_AHM7KMBGXK/HM7KMBGXK"
outDir="/data/projects/220311_Brandt_Berres_MiedIII_scCITEseq/analysis/CITE_count"
baseDir="/data/projects/220311_Brandt_Berres_MiedIII_scCITEseq/analysis"
refGenome="/data/shared_env/10xGenomics/refdata-gex-GRCh38-2020-A"
refFeature="/data/projects/220311_Brandt_Berres_MiedIII_scCITEseq/analysis/EB220311_feature_ref.csv"

# echo "$refFeature"
# echo "$refGenome"


mkdir -p $outDir



cd $fastqDir

samples=`find . -type f  -exec basename "{}" \; | cut -d '_' -f1 | sort -u`

cd $outDir

for sampleid in $samples
do
			LibCsv="$baseDir"/"$sampleid"_config.csv
#			echo "$LibCsv"
#			echo "$sampleid"

			cellranger count --id=$sampleid \
					--libraries=$LibCsv \
					--feature-ref=$refFeature \
					--transcriptome=$refGenome \
					--localcores=15 \
					--localmem=120 \

done
