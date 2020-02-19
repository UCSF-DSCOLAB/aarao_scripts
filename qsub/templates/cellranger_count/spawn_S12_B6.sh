#!/bin/bash

set -e
set -o nounset
export MEMORY=384
export TRANSCRIPTOME=/data/shared/krummellab/ipi/data/refs/10x/HetGla02

base_dir=/data/shared/krummellab/immunox/XHM/sequencing/fastqs
log_dir=/data/shared/krummellab/arrao/logs
mkdir -p ${log_dir}

tree -fid -L 3 ${base_dir} | grep fastqs$ | grep "B6\|S[12]" | while read fastq_dir
do
    export OUTDIR=`echo $fastq_dir | sed s/"fastqs"/"processed"/g`
    export FASTQDIR=${fastq_dir}
    export SAMPLE=`echo $(basename ${fastq_dir}) | sed s/"_fastqs"//g`
    export SAMPLEID=${SAMPLE}
    echo $OUTDIR $FASTQDIR $SAMPLE $SAMPLEID
    qsub -V \
         -e ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N cellranger_count_${SAMPLEID} \
         -l nodes=1:ppn=64 \
         -l vmem=${MEMORY}gb \
         run.sh
done
