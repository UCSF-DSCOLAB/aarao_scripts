#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/arrao/MarkDuplicates_${SAMPLE} && cd /scratch/arrao/MarkDuplicates_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/MarkDuplicates_${SAMPLE} ; }" EXIT

out_dir=$(dirname ${BAMFILE})
out_base=$(basename ${BAMFILE%.bam})_DupMarked
java -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp \
    -jar /data/shared/krummellab/ipi/software/picard-2.18.14/picard.jar \
    MarkDuplicates \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/scratch/arrao/${SAMPLE}_picardtemp \
    REMOVE_DUPLICATES=false \
    AS=true \
    I=${BAMFILE} \
    O=${out_base}.bam \
    M=${out_base}.duplication_metrics

mv /scratch/arrao/MarkDuplicates_${SAMPLE}/${out_base}* ${out_dir}/