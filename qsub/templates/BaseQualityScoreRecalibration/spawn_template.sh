#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=1
#@ MEMREQS vmem=150gb
#@ REFERENCEFASTA /data/shared/krummellab/ipi/data/refs/hg38_files/hg38.fa
#@ DBSNP /data/shared/krummellab/ipi/data/databases/dbsnp/hg38/b151/00-All.vcf

set -e
set -o nounset

export _memory_=`echo "$(echo "<_#MEMREQS#_>" | sed 's/[^0-9]*//g')*0.9 / 1" | bc`
export _reference_fa_=<_#REFERENCEFASTA#_>
export _dbsnp_=<_#DBSNP#_>

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
  
  _bamfile_="${bam_dir}/$(ls ${bam_dir} | grep "bam$")"
  export _sample_=${line}

  qsub -v MEMORY=${_memory_},REFERENCE_FA=${_reference_fa_},DBSNP=${_dbsnp_},BAMFILE=${_bamfile_},SAMPLE=${_sample_} \
       -e ${log_dir}/BQSR_${_sample_}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/BQSR_${_sample_}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N BQSR_${_sample_} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done  


