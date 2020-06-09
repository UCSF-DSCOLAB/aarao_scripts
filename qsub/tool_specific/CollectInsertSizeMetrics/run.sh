#!/bin/bash

# For Rscript histogram plotting. Honestly doesn't really matter what env is used
# as long as it has R and ggplot2.
source /krummellab/data1/ipi/software/miniconda3/etc/profile.d/conda.sh
conda activate r_seurat_3_6

set -e
set -o nounset

module load CBC gatk/${GATKVERSION}   

mkdir /scratch/${USER}/ISizeMetrics_${SAMPLE} && cd /scratch/${USER}/ISizeMetrics_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/ISizeMetrics_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT


out_dir=$(dirname ${SAMFILE})
gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   CollectInsertSizeMetrics \
   -I ${SAMFILE} \
   -O insert_size_metrics.txt \
   -H insert_size_histogram.pdf

mv /scratch/${USER}/ISizeMetrics_${SAMPLE}/insert_size_metrics* ${out_dir}/
