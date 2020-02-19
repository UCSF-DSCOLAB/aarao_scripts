#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=24
#@ MEMREQS vmem=150gb
#@ MINREADLEN 30
#@ KIT nextera

set -e
set -o nounset

base_dir=<__BASEDIR__>
log_dir=<_#LOGDIR#_>
out_dir=<__OUTDIR__>
kit=<_#KIT#_>
mkdir -p ${log_dir} ${out_dir}

export OUTDIR=$out_dir

if [ "$kit" == "truseq" ]
then
    export READ1ADAPTER=AGATCGGAAGAG
    export READ2ADAPTER=AGATCGGAAGAG
    export NEXTSEQTRIM=''
elif [ "$kit" == "nextera" ]
then
    export READ1ADAPTER=CTGTCTCTTATACACATCT
    export READ2ADAPTER=CTGTCTCTTATACACATCT
    export NEXTSEQTRIM="--nextseq-trim=20"
else
    echo "Kit can only be nextera or truseq"
    exit 1
fi

ls ${base_dir} | while read line
do
  if [ ! -d ${base_dir}/${line} ]
  then
    continue
  fi
  export FQDIR=${base_dir}/${line}
  if [ $(ls ${FQDIR} | grep -c "fastq\|fq") -ne 2 ]
  then
    echo "Cannot process ${FQDIR} since it does not have 2 fastq files."
    continue
  fi
  fq_files=($(ls ${FQDIR} | sort | grep "fastq\|fq"))
  f_diff=($(cmp -bl <(echo ${fq_files[0]}) <(echo  ${fq_files[1]}) || :))
  if [ ${#f_diff[@]} -ne 5 ]
  then
    echo "Couldn't decompose fastqs in ${FQDIR} into 'PREFIX[12]SUFFIX' because there are > 1 differences between the 2 fastq names."
    continue
  elif [ ${f_diff[2]}${f_diff[4]} -ne 12 ]
  then
    echo "Couldn't decompose fastqs in ${FQDIR} into 'PREFIX[12]SUFFIX' because the difference between the names was not 1/2."
    continue
  fi
  export FQPREFIX=$(echo ${fq_files[0]} | cut -c -$((${f_diff[0]} - 1)))
  export FQSUFFIX=$(echo ${fq_files[0]} | cut -c $((${f_diff[0]} + 1))-)
  export MINREADLEN=<_#MINREADLEN#_>
  export SAMPLE=${line}

  qsub -V \
       -e ${log_dir}/cutadapt_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       -o ${log_dir}/cutadapt_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       -N cutadapt_${SAMPLE} \
       -l <_#NODEREQS#_> \
       -l <_#MEMREQS#_> \
       run.sh
done


