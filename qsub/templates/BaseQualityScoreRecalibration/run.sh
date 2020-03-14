#!/bin/bash
set -e
set -o nounset

module load CBC gatk/${GATKVERSION}

mkdir /scratch/arrao/BQSR_${SAMPLE} && cd /scratch/arrao/BQSR_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/BQSR_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_BQSRd

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    BaseRecalibrator\
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --known-sites ${DBSNP} \
    --output recal_data.table

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    ApplyBQSR \
    --create-output-bam-index true \
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --bqsr-recal-file recal_data.table \
    --output ${out_base}.bam

mv /scratch/arrao/BQSR_${SAMPLE}/recal_data.table ${out_dir}/
mv /scratch/arrao/BQSR_${SAMPLE}/${out_base}* ${out_dir}/
