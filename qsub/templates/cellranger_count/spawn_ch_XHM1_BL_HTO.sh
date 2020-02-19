#!/bin/bash

set -e
set -o nounset

MEMORY=vmem=250gb
MEMORY=${MEMORY#vmem=}
MEMORY=${MEMORY%gb}

libraries_csv=/krummellab/data1/immunox/XHM/sequencing/fastqs/pooled/BL_10X_libraries.csv
log_dir=/krummellab/data1/arrao/logs
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

qsub -v "MEMORY=${MEMORY},FEATURE_REF=/krummellab/data1/ipi/data/refs/10x/biolegend_totalseq_hashtags.csv,CHEMISTRY=auto,TRANSCRIPTOME=/krummellab/data1/ipi/data/refs/10x/GRCm38,OUTDIR=${OUTDIR},LIBRARIES_CSV=${libraries_csv},SAMPLEID=${sample//./_}" \
     -e ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/cellranger_count_${sample}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_count_${sample} \
     -l nodes=1:ppn=32 \
     -l vmem=250gb \
     run_ch.sh
