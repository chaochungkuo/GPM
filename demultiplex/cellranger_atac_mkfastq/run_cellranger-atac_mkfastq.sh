# Please make sure cellranger-atac-2.1.0 is available in your environment
# You should review the paramenters and run:
# screen -S cellranger
# bash run_cellranger-atac_mkfastq.sh

# Please execute this command in the directory OUTPUT_DIR
PATH_CELLRANGER_ATAC mkfastq --id=mkfastq --localcores=N_CORES \
                   --run=BCL_PATH \
                   --csv=./samplesheet.csv

bash run_merge_lanes.sh