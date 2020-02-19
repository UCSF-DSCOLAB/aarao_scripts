#!/bin/bash
set -e
set -o nounset

fq_dir=/data/shared/krummellab/ipi/fastqs/Plate27Rnaseq
out_dir=/data/shared/krummellab/ipi/qc
tinfoil_dir=${out_dir}/tinfoil/Plate27Rnaseq
fastqc_dir=${out_dir}/fastqc/Plate27Rnaseq
log_dir=${out_dir}/logs/Plate27Rnaseq
mkdir -p ${log_dir} ${tinfoil_dir} ${fastqc_dir}

tree -fi ${fq_dir} | grep fastq.gz | sed s/".fastq.gz"//g | while read line
do
  export FASTQ_FILE=${line}.fastq.gz
  export SAMPLE=$(basename ${line})
  export TINFOILDIR=${tinfoil_dir}
  export FASTQCDIR=${fastqc_dir}
  export READLEN=151
  qsub -V \
       -e ${log_dir}/tinfoil_${SAMPLE}.err \
       -o ${log_dir}/tinfoil_${SAMPLE}.out \
       -N qc_27_${SAMPLE} \
       -l nodes=1:ppn=1 \
       -l vmem=50gb \
       run.sh
done

