#!/bin/bash
set -e
set -o nounset

mkdir /scratch/arrao/bwa_${SAMPLE} && cd /scratch/arrao/bwa_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/bwa_${SAMPLE} ; }" EXIT

mkdir -p ${OUTDIR}

bwa mem -t ${PBS_NUM_PPN} -v 1 /data/shared/krummellab/ipi/data/refs/bwa/hg38 ${FQDIR}/${FQPREFIX}1${FQSUFFIX} ${FQDIR}/${FQPREFIX}2${FQSUFFIX} > ${SAMPLE}.sam
samtools view -bS -o ${SAMPLE}.bam ${SAMPLE}.sam && rm ${SAMPLE}.sam
samtools sort ${SAMPLE}.bam ${SAMPLE}_sorted
samtools index ${SAMPLE}_sorted.bam
mv ${SAMPLE}_sorted* ${OUTDIR}/ 
