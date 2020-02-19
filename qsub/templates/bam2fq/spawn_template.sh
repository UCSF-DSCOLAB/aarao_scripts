#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ OUTDIR <__BAMDIR__>/fastqs
#@ NODEREQS nodes=1:ppn=1
#@ MEMREQS vmem=150gb

set -e
set -o nounset

bam_dir=<__BAMDIR__>
out_dir=<_#OUTDIR#_>
log_dir=<_#LOGDIR#_>
mkdir -p ${log_dir} ${out_dir}

tree -fi ${bam_dir} | grep bam$ | sed s/".bam"//g | while read line
do
  export BAMFILE=${line}.bam
  export SAMPLE=$(basename ${line})
  export OUTDIR=${out_dir}
  qsub -V \
       -e ${log_dir}/bam2fq_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/bam2fq_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N bam2fq_${SAMPLE} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done