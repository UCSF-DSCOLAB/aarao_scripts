#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/${USER}/MarkDuplicates_${SAMPLE}
cd /scratch/${USER}/MarkDuplicates_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/MarkDuplicates_${SAMPLE} ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_DupMarked

cp ${SAMFILE} .
java -XX:ParallelGCThreads=${PBS_NUM_PPN} \
    -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp \
    -jar /krummellab/data1/ipi/software/picard/${PICARDVERSION}/picard.jar \
    MarkDuplicates \
    VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY} \
    TMP_DIR=/scratch/${USER}/${SAMPLE}_picardtemp \
    REMOVE_DUPLICATES=${REMOVE_DUPLICATES} \
    REFERENCE_SEQUENCE=${GENOMEREF} \
    ASSUME_SORTED=${ASSUME_SORTED} \
    CREATE_INDEX=${CREATE_INDEX} \
    INPUT=$(basename ${SAMFILE}) \
    OUTPUT=${out_base}.bam \
    METRICS_FILE=${out_base}.duplication_metrics

mv /scratch/${USER}/MarkDuplicates_${SAMPLE}/${out_base}* ${out_dir}/
