#!/bin/bash
set -e
set -o nounset

module load CBC gatk/${GATKVERSION}

mkdir /scratch/${USER}/BQSR_${SAMPLE} && cd /scratch/${USER}/BQSR_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/BQSR_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_BQSRd

gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    BaseRecalibrator\
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --known-sites ${DBSNP} \
    --output ${out_base}_recal_data.table

gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    ApplyBQSR \
    --create-output-bam-index true \
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --bqsr-recal-file ${out_base}_recal_data.table \
    --output ${out_base}.bam

mv /scratch/${USER}/BQSR_${SAMPLE}/${out_base}* ${out_dir}/
