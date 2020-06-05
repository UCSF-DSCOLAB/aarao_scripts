#!/bin/bash
set -e
set -o nounset

module load CBC gatk/${GATKVERSION}

mkdir /scratch/arrao/IndelRealignment_${SAMPLE} && cd /scratch/arrao/IndelRealignment_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/IndelRealignment_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_IndelRealigned

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    LeftAlignIndels \
    --create-output-bam-index true \
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --OUTPUT ${out_base}.bam

mv /scratch/arrao/IndelRealignment_${SAMPLE}/${out_base}* ${out_dir}/
