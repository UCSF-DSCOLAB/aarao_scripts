#!/bin/bash
#@ LOGDIR /krummellab/data1/arrao/logs
#@ NODEREQS nodes=1:ppn=64
#@ MEMREQS vmem=500gb

set -e
set -o nounset

source ~/.bash_functions

NODEREQS="<_#NODEREQS#_>"
MEMREQS="<_#MEMREQS#_>"
MEMORY=${MEMREQS#vmem=}
MEMORY=${MEMORY%gb}

sample_sheet=<__SAMPLESHEET__>
bcl_dir=<__BCLDIR__>
# Strip trailing slash
bcl_dir=${bcl_dir%/}

log_dir=<_#LOGDIR#_>
# Strip trailing slash
log_dir=${log_dir%/}
mkdir -p ${log_dir}

OUTDIR=${bcl_dir}_fastqs
uuid=`randomstr 10`
echo qsub -v "UUID=${uuid},MEMORY=${MEMORY},OUTDIR=${OUTDIR},SAMPLESHEET=${sample_sheet},BCLDIR=${bcl_dir}" \
     -e ${log_dir}/cellranger_vdj_${uuid}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/cellranger_vdj_${uuid}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_vdj_${uuid} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     -l feature=n34 \
     run.sh