#!/bin/bash

source /krummellab/data1/ipi/software/miniconda3/etc/profile.d/conda.sh
conda activate r_seurat_3_6

set -e
set -o nounset

if [ "${WORKING_FOLDER}" == "TMPDIR" ]
then
  WORKING_FOLDER_STRING="WORKING_FOLDER=${TMPDIR}"
else
  WORKING_FOLDER_STRING="WORKING_FOLDER=${WORKING_FOLDER}"
fi

Rscript /krummellab/data1/${USER}/aarao_scripts/R/process_10x_with_Seurat.R SAMPLE_YML=${SAMPLE_YML} RSCRIPTS_DIR=${RSCRIPTS_DIR} GENESET_DIR=${GENESET_DIR} ${WORKING_FOLDER_STRING}
