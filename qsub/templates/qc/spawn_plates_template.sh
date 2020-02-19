#!/bin/bash
#@ READLEN 101
set -e
set -o nounset

fq_dir=/cbc2/data2/samadb/IPI/fastqs/Plate<__PLATE__>Rnaseq
out_dir=/data/shared/krummellab/ipi/qc/Plate<__PLATE__>Rnaseq
log_dir=${out_dir}/logs/Plate<__PLATE__>Rnaseq
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
       -N qc_<__PLATE__>_${SAMPLE} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done
