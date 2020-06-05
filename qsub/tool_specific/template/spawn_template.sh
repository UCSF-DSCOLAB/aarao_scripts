#!/bin/bash
#@ LOGDIR /data/shared/krummellab/arrao/logs
#@ NODEREQS nodes=1:ppn=1
#@ MEMREQS vmem=150gb
set -e
set -o nounset

qsub -V \
     -e ${log_dir}/xxx_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${log_dir}/xxx_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N xxx_${SAMPLE} \
     -l <_#NODEREQS#_> \
     -l <_#MEMREQS#_> \
     run.sh
