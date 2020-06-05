#!/bin/bash
#@ LOGDIR /krummellab/data1/arrao/logs
#@ NODEREQS nodes=1:ppn=64
#@ MEMREQS vmem=500gb
#@ TRANSCRIPTOME /krummellab/data1/ipi/data/refs/10x/hg38
#@ CHEMISTRY auto
#@ FEATUREREF /krummellab/data1/ipi/data/refs/10x/biolegend_totalseq_hashtags.csv

set -e
set -o nounset

MEMORY=<_#MEMREQS#_>
MEMORY=${MEMORY#vmem=}
MEMORY=${MEMORY%gb}

libraries_csv=<__LIBRARIESCSV__>
log_dir=<_#LOGDIR#_>
# Strip trailing slash
log_dir=${log_dir%/}
mkdir -p ${log_dir}

if [ ! -f ${libraries_csv} ]
then
    echo "Invalid libraries_csv. Does not exist."
    exit 1
fi

gex_dir=`grep "Gene Expression" ${libraries_csv} | awk -F"," '{print $1}'`
sample=`grep "Gene Expression" ${libraries_csv} | awk -F"," '{print $2}'`

if [[ $gex_dir =~ "/fastqs/" ]]
then
    OUTDIR=`echo ${gex_dir} | sed s/"\/fastqs\/"/"\/processed\/"/g`
else
    OUTDIR=${gex_dir}/processed
fi
mkdir -p ${OUTDIR}

qsub -v "MEMORY=${MEMORY},FEATURE_REF=<_#FEATUREREF#_>,CHEMISTRY=<_#CHEMISTRY#_>,TRANSCRIPTOME=<_#TRANSCRIPTOME#_>,OUTDIR=${OUTDIR},LIBRARIES_CSV=${libraries_csv},SAMPLEID=${sample//./_}" \
     -e ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_count_${sample} \
     -l <_#NODEREQS#_> \
     -l <_#MEMREQS#_> \
     run_ch.sh
