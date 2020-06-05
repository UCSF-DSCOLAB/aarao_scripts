#!/bin/bash
#@ CHEM V2
#@ CELLS 50000
#@ CSTAGS /data/shared/krummellab/ipi/data/refs/10x/biolegend_totalseq_citeseqcount_tags.csv
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=10
#@ MEMREQS vmem=150gb
set -e
set -o nounset

base_dir=<__BASEDIR__>
log_dir=<_#LOGDIR#_>
mkdir -p ${log_dir}

ls ${base_dir} | grep _R1_.*.fastq.gz | while read fq1_file
do
    prefix=`echo ${fq1_file} | sed s/"_R1_.*"//g`
    sample=<__SAMPLE__>_${prefix}
    fq_files=($(ls ${base_dir}/${prefix}* | grep -v '_I1_'))
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
    fqprefix=$(echo ${fq_files[0]} | cut -c -$((${f_diff[0]} - 1)))
    fqsuffix=$(echo ${fq_files[0]} | cut -c $((${f_diff[0]} + 1))-)
    fq1=${fqprefix}1${fqsuffix}
    fq2=${fqprefix}2${fqsuffix}

    qsub -v "FQ1=${fq1},FQ2=${fq2},CHEM=<_#CHEM#_>,CELLS=<_#CELLS#_>,OUTDIR=<__OUTDIR__>,SAMPLE=${sample},CSTAGS=<_#CSTAGS#_>" \
         -e ${log_dir}/citeseq_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${log_dir}/citeseq_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N citeseq_count_${sample} \
         -l <_#NODEREQS#_> \
         -l <_#MEMREQS#_> \
         run.sh
done
