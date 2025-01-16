#!/bin/bash

PATH_NEXTFLOW run nf-core/demultiplex -r 1.5.4 \
   -profile docker -c nextflow.config \
   --input samplesheet.csv \
   --outdir results