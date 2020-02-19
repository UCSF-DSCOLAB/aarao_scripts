#!/bin/bash

set -e
set -o nounset

source ~/.bash_functions

NODEREQS="nodes=1:ppn=32"
MEMREQS="vmem=250gb"
MEMORY=${MEMREQS#vmem=}
MEMORY=${MEMORY%gb}

base_dir=/krummellab/data1/ipi/sequencing/fastqs/10X/5prime_TCR/CRC
# Strip trailing slash
base_dir=${base_dir%/}
log_dir=/krummellab/data1/arrao/logs
# Strip trailing slash
log_dir=${log_dir%/}
mkdir -p ${log_dir}

ls ${base_dir} | while read sample
do
    fastq_dir=${base_dir}/${sample}
    if [ ! -d ${fastq_dir} ]
    then
        continue
    fi
    OUTDIR=`echo ${fastq_dir} | sed s/"\/fastqs\/"/"\/processed\/"/g`

    uuid=`randomstr 10`
    qsub -v "UUID=${uuid},MEMORY=${MEMORY},REFERENCE=/krummellab/data1/ipi/data/refs/10x/GRCh38_VDJ,OUTDIR=${OUTDIR},FASTQDIR=${fastq_dir},SAMPLEID=${sample//./_},CHAIN=TR" \
         -e ${log_dir}/cellranger_vdj_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${log_dir}/cellranger_vdj_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N cellranger_vdj_${sample} \
         -l ${NODEREQS} \
         -l ${MEMREQS} \
         run.sh
done
