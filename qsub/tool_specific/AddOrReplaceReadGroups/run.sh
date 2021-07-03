#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/${USER}/AddRepRGs_${SAMPLE}
cd /scratch/${USER}/AddRepRGs_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/AddRepRGs_${SAMPLE} ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_AddRepRGd
java -XX:ParallelGCThreads=${PBS_NUM_PPN} \
    -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp \
    -jar /krummellab/data1/ipi/software/picard/${PICARDVERSION}/picard.jar \
    AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY} \
    TMP_DIR=/scratch/${USER}/${SAMPLE}_picardtemp \
    ASSUME_SORTED=${ASSUME_SORTED} \
    CREATE_INDEX=${CREATE_INDEX} \
    INPUT=${SAMFILE} \
    OUTPUT=${out_base}.bam \
    SORT_ORDER=coordinate \
    RGID=${RGID} \
    RGLB=${RGLB} \
    RGPL=${RGPL} \
    RGPU=${RGPU} \
    RGSM=${RGSM}   

mv /scratch/${USER}/AddRepRGs_${SAMPLE}/${out_base}* ${out_dir}/
