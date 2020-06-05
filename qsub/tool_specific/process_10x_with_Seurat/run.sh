#!/bin/bash

source /home/arrao/miniconda3/etc/profile.d/conda.sh
conda activate r_monocle2

set -e
set -o nounset

# None of thiw works right now!
#start=`grep -n "args = list" /krummellab/data1/arrao/scripts/R/process_10x_with_Seurat.R | awk -F ":" '{print $1}'`
#end=`grep -n "argsClasses = " /krummellab/data1/arrao/scripts/R/process_10x_with_Seurat.R | awk -F ":" '{print $1}'`
#expected_args=(`tail -n+$((${start} + 1)) /krummellab/data1/arrao/scripts/R/process_10x_with_Seurat.R | 
#               head -$((${end} - ${start} - 1 )) | 
#               grep "=" | 
#               sed s/"=.*"//g`)

#argstring=""
#for s in ${expected_args[@]}
#do
#    argstring+="${s}=`printenv ${s}` "
#done

Rscript /krummellab/data1/arrao/scripts/R/process_10x_with_Seurat.R SAMPLE_YML=${SAMPLE_YML} OUT_FOLDER=${OUT_FOLDER}