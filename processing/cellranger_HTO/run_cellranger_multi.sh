#!/bin/bash

Local_Cores=50

cellranger="/data/shared/10xGenomics/cellranger-9.0.0/bin/cellranger"

${cellranger} multi --id="pool1" --csv="multi_.csv" --localcores=$Local_Cores
