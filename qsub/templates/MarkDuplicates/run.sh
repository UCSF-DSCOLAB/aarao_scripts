#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/arrao/MarkDuplicates_${SAMPLE} && cd /scratch/arrao/MarkDuplicates_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/MarkDuplicates_${SAMPLE} ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_DupMarked
java -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp \
    -jar /krummellab/data1/ipi/software/picard/${PICARDVERSION}/picard.jar \
    MarkDuplicates \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/scratch/arrao/${SAMPLE}_picardtemp \
    REMOVE_DUPLICATES=false \
    REFERENCE_SEQUENCE=${GENOMEREF} \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    INPUT=${SAMFILE} \
    OUTPUT=${out_base}.bam \
    METRICS_FILE=${out_base}.duplication_metrics

mv /scratch/arrao/MarkDuplicates_${SAMPLE}/${out_base}* ${out_dir}/