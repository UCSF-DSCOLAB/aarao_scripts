#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/arrao/AddRepRGs_${SAMPLE} && cd /scratch/arrao/AddRepRGs_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/AddRepRGs_${SAMPLE} ; }" EXIT

out_dir=$(dirname ${BAMFILE})
out_base=$(basename ${BAMFILE%.bam})_RGFixed
java -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp \
    -jar /data/shared/krummellab/ipi/software/picard-2.18.14/picard.jar \
    AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true \
    I=${BAMFILE} \
    O=${out_base}.bam \
    SO=coordinate \
    ID=1 \
    LB=${SAMPLE} \
    PL=ILLUMINA \
    PU=12345 \
    SM=${SAMPLE}_exome 

mv /scratch/arrao/AddRepRGs_${SAMPLE}/${out_base}* ${out_dir}/