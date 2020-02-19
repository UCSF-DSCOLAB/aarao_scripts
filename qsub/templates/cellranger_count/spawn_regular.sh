#!/bin/bash

set -e
set -o nounset
export MEMORY=190
export TRANSCRIPTOME=/data/shared/krummellab/ipi/data/refs/10x/hg38

base_dir=/data/shared/krummellab/immunox/XNI
log_dir=/data/shared/krummellab/arrao/logs
mkdir -p ${log_dir}

tree -fid -L 2 ${base_dir} | grep fastqs$ | while read fastq_dir
do
    export OUTDIR=$(dirname ${fastq_dir})/results
    ls ${fastq_dir} | grep BRN_6_S2 | while read sample
    do
        export FASTQDIR=${fastq_dir}/${sample}
        export SAMPLE=`sed s/"_fastqs"/""/g <(echo ${sample})`
        export SAMPLEID=${SAMPLE}_default_cells
        qsub -V \
             -e ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
             -o ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
             -N cellranger_count_${SAMPLEID} \
             -l nodes=1:ppn=24 \
             -l vmem=190gb \
             run.sh
    done
done
