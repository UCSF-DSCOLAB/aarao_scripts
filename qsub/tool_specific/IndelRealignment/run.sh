#!/bin/bash
set -e
set -o nounset

module load CBC gatk/${GATKVERSION}

mkdir /scratch/${USER}/IndelRealignment_${SAMPLE}
cd /scratch/${USER}/IndelRealignment_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/IndelRealignment_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_IndelRealigned

gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    LeftAlignIndels \
    --create-output-bam-index true \
    --reference ${GENOMEREF} \
    --input ${SAMFILE} \
    --OUTPUT ${out_base}.bam

mv /scratch/${USER}/IndelRealignment_${SAMPLE}/${out_base}* ${out_dir}/
