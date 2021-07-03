#!/bin/bash

set -e
set -o nounset

module load CBC gatk/${GATKVERSION}   

mkdir /scratch/${USER}/CollectHsMetrics_${SAMPLE}
cd /scratch/${USER}/CollectHsMetrics_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/CollectHsMetrics_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp /scratch/${USER}/${SAMPLE}_picardtemp ; }" EXIT


extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
java -XX:ParallelGCThreads=${PBS_NUM_PPN} \
    -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp \
    -jar /krummellab/data1/ipi/software/picard/${PICARDVERSION}/picard.jar \
    CollectHsMetrics \
    VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY} \
    TMP_DIR=/scratch/${USER}/${SAMPLE}_picardtemp \
    INPUT=${SAMFILE} \
    OUTPUT=$(basename ${SAMFILE%.${extension}})_hs_metrics.txt \
    REFERENCE_SEQUENCE=${GENOMEREF} \
    BAIT_INTERVALS=${BAIT_INTERVALS} \
    TARGET_INTERVALS=${TARGET_INTERVALS}

mv /scratch/${USER}/CollectHsMetrics_${SAMPLE}/*_hs_metrics.txt ${out_dir}/
