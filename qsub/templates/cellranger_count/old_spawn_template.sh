#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=24
#@ MEMREQS vmem=150gb
#@ MEMORY 150
#@ TRANSCRIPTOME /data/shared/krummellab/ipi/data/refs/10x/hg38

set -e
set -o nounset
export MEMORY=<_#MEMORY#_>
export TRANSCRIPTOME=<_#TRANSCRIPTOME#_>

base_dir=<__BASEDIR__>
log_dir=<_#LOGDIR#_>
mkdir -p ${log_dir}

tree -fid -L 2 ${base_dir} | grep fastqs$ | while read fastq_dir
do
    export OUTDIR=$(dirname ${fastq_dir})/results
    ls ${fastq_dir} | while read sample
    do
        export FASTQDIR=${fastq_dir}/${sample}
        export SAMPLE=${sample}
        export SAMPLEID=${SAMPLE}_default_cells
        qsub -V \
             -e ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
             -o ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
             -N cellranger_count_${SAMPLEID} \
             -l <_#NODEREQS#_> \
             -l <_#MEMREQS#_> \
             run.sh
    done
done
