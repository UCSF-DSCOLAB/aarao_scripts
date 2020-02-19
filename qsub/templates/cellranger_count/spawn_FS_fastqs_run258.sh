#!/bin/bash

set -e
set -o nounset

MEMORY=vmem=250gb
MEMORY=${MEMORY#vmem=}
MEMORY=${MEMORY%gb}

fastq_dir=/krummellab/data1/immunox/XCAN1/experiments/fastqs/10x_3prime/FS_fastqs_run258
# Strip trailing slash
fastq_dir=${fastq_dir%/}
log_dir=/krummellab/data1/arrao/logs
# Strip trailing slash
log_dir=${log_dir%/}
mkdir -p ${log_dir}

sample=$(basename ${fastq_dir})
if [ ! -d ${fastq_dir} ]
then
    echo "Invalid fastq dir. Not a directory."
    exit 1
fi

if [[ $fastq_dir =~ "/fastqs/" ]]
then
    OUTDIR=`echo ${fastq_dir} | sed s/"\/fastqs\/"/"\/processed\/"/g`
else
    OUTDIR=${fastq_dir}/processed
fi
mkdir -p ${OUTDIR}

qsub -v "MEMORY=${MEMORY},CHEMISTRY=auto,TRANSCRIPTOME=/krummellab/data1/ipi/data/refs/10x/hg38,OUTDIR=${OUTDIR},FASTQDIR=${fastq_dir},SAMPLEID=${sample//./_}" \
     -e ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_count_${sample} \
     -l nodes=1:ppn=32 \
     -l vmem=250gb \
     run.sh

