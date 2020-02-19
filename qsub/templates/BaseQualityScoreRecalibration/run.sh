#!/bin/bash
set -e
set -o nounset

module load CBC gatk/4.0.2.1

mkdir /scratch/arrao/BQSR_${SAMPLE} && cd /scratch/arrao/BQSR_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/BQSR_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT

out_dir=$(dirname ${BAMFILE})
out_base=$(basename ${BAMFILE%.bam})_BQSRd

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    BaseRecalibrator\
    --reference ${REFERENCE_FA} \
    --input ${BAMFILE} \
    --known-sites ${DBSNP} \
    --output recal_data.table

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    ApplyBQSR \
    --reference ${REFERENCE_FA} \
    --input ${BAMFILE} \
    --bqsr-recal-file recal_data.table \
    --output ${out_base}.bam

mv /scratch/arrao/BQSR_${SAMPLE}/recal_data.table ${out_dir}/
mv /scratch/arrao/BQSR_${SAMPLE}/${out_base}* ${out_dir}/
