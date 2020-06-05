#!/bin/bash
#@ READLEN 101
set -e
set -o nounset

fq_dir=<__FASTQDIR__>
out_dir=<__FASTQDIR__>/qc
log_dir=<_#LOGDIR#_>
mkdir -p ${log_dir} ${out_dir}

tree -fi ${fq_dir} | grep fastq.gz | sed s/".fastq.gz"//g | while read line
do
  export FASTQ_FILE=${line}.fastq.gz
  export SAMPLE=$(basename ${line})
  export OUTDIR=${out_dir}
  export READLEN=<_#READLEN#_>
  qsub -V \
       -e ${log_dir}/qc_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/qc_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N qc_${SAMPLE} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done
