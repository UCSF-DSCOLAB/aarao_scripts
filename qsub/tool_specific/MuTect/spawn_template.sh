#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=1
#@ MEMREQS vmem=150gb

set -e
set -o nounset

memory=`echo "$(echo "<_#MEMREQS#_>" | sed 's/[^0-9]*//g')*0.90 / 1" | bc`

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

  qsub -v "MEMORY=${memory},BAMFILE=${bamfile},SAMPLE=${sample}" \
       -e ${log_dir}/MuTect_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/MuTect_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N MuTect_${sample} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done  


