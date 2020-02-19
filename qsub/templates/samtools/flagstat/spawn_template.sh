#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=1
#@ MEMREQS vmem=150gb

set -e
set -o nounset

base_dir=<__BASEDIR__>
log_dir=<_#LOGDIR#_>
mkdir -p ${log_dir}

ls ${base_dir} | while read line
do
  if [ ! -d ${base_dir}/${line} ]
  then
    continue
  fi
  bam_dir=${base_dir}/${line}
  if [ $(ls ${bam_dir} | grep -c "bam$") -ne 1 ]
  then
    echo "Cannot process ${bam_dir} since it does not have exactly 1 bam."
    continue
  fi
  
  bamfile=${bam_dir}/$(ls ${bam_dir} | grep "bam$")
  sample=${line}

  qsub -v "BAMFILE=${bamfile}" \
       -e ${log_dir}/samtools_flagstat_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/samtools_flagstat_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N samtools_flagstat_${sample} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done  


