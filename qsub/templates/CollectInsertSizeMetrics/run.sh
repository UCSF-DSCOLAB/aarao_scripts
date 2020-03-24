#!/bin/bash

# For Rscript histogram plotting
source /home/arrao/miniconda3/etc/profile.d/conda.sh
conda activate r_monocle2

set -e
set -o nounset

module load CBC gatk/${GATKVERSION}   

mkdir /scratch/arrao/ISizeMetrics_${SAMPLE} && cd /scratch/arrao/ISizeMetrics_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/ISizeMetrics_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT


out_dir=$(dirname ${SAMFILE})
gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   CollectInsertSizeMetrics \
   -I ${SAMFILE} \
   -O insert_size_metrics.txt \
   -H insert_size_histogram.pdf

mv /scratch/arrao/ISizeMetrics_${SAMPLE}/insert_size_metrics* ${out_dir}/
