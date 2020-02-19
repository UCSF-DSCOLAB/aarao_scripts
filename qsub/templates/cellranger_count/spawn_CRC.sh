#!/bin/bash

set -e
set -o nounset

MEMORY=vmem=250gb
MEMORY=${MEMORY#vmem=}
MEMORY=${MEMORY%gb}

base_dir=/krummellab/data1/ipi/sequencing/fastqs/10X/5prime_GEX/CRC
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
    export OUTDIR=`echo ${fastq_dir} | sed s/"\/fastqs\/"/"\/processed\/"/g`

    echo qsub -v "MEMORY=${MEMORY},TRANSCRIPTOME=/krummellab/data1/ipi/data/refs/10x/hg38,OUTDIR=${OUTDIR},FASTQDIR=${fastq_dir},SAMPLEID=${sample//./_}" \
         -e ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N cellranger_count_${sample} \
         -l nodes=1:ppn=32 \
         -l vmem=250gb \
         run.sh
done
