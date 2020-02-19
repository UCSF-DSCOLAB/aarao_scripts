#!/bin/bash
#@ LOGDIR /krummellab/data1/arrao/logs
#@ NODEREQS nodes=1:ppn=32
#@ MEMREQS vmem=250gb
#@ TRANSCRIPTOME /krummellab/data1/ipi/data/refs/10x/hg38
#@ CHEMISTRY auto

set -e
set -o nounset

MEMORY=<_#MEMREQS#_>
MEMORY=${MEMORY#vmem=}
MEMORY=${MEMORY%gb}

fastq_dir=<__FASTQDIR__>
# Strip trailing slash
fastq_dir=${fastq_dir%/}
log_dir=<_#LOGDIR#_>
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

qsub -v "MEMORY=${MEMORY},CHEMISTRY=<_#CHEMISTRY#_>,TRANSCRIPTOME=<_#TRANSCRIPTOME#_>,OUTDIR=${OUTDIR},FASTQDIR=${fastq_dir},SAMPLEID=${sample//./_}" \
     -e ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_count_${sample} \
     -l <_#NODEREQS#_> \
     -l <_#MEMREQS#_> \
     run.sh

