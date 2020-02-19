#!/bin/bash

set -e
set -o nounset
export MEMORY=384
export TRANSCRIPTOME=/data/shared/krummellab/ipi/data/refs/10x/hg38

base_dir=/cbc/samadb/ipi_fastqs/10X_3prime/hashed/test1
log_dir=/krummellab/data1/arrao/logs
mkdir -p ${log_dir}

tree -fid -L 1 ${base_dir} | grep fastqs$ | while read fastq_dir
do
    export FASTQDIR=${fastq_dir}
    export SAMPLE=`echo $(basename ${fastq_dir}) | sed s/"_fastqs"//g`
    export SAMPLEID=${SAMPLE}
    export OUTDIR=${base_dir}/output/${SAMPLE}
    echo $OUTDIR $FASTQDIR $SAMPLE $SAMPLEID
    qsub -V \
         -e ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${log_dir}/cellranger_count_${SAMPLEID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N cellranger_count_${SAMPLEID} \
         -l nodes=1:ppn=64 \
         -l vmem=${MEMORY}gb \
         run_ch.sh
done
