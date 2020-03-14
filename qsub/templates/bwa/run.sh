#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/ipi/ipi_usr/SOURCE_THIS 

mkdir /scratch/arrao/bwa_${SAMPLE} && cd /scratch/arrao/bwa_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/bwa_${SAMPLE} ; }" EXIT

getRG(){
 zcat ${1} | 
 head -1 | 
 awk -F":" -v SAMPLE=${2} '{print "@RG\\tID:"$3"."$2"."SAMPLE"."$4"\\tPL:ILLUMINA\\tSM:"SAMPLE"\\tPU:"$3"."$2"."SAMPLE"."$4"\\tBC:"$NF}'
}

RG=`getRG ${FQ1} ${SAMPLE}`

bwa mem -t ${PBS_NUM_PPN} \
        -v 1 \
        -R ${RG} \
        ${BWAREF} \
        ${FQ1} \
        ${FQ2} > ${SAMPLE}.sam

samtools view -b \
              -o ${SAMPLE}.bam \
              -@ ${PBS_NUM_PPN} \
              ${SAMPLE}.sam \
              && rm ${SAMPLE}.sam

samtools sort ${SAMPLE}.bam \
              -@ ${PBS_NUM_PPN} \
              -o ${SAMPLE}_sorted.bam \
              && rm ${SAMPLE}.bam

samtools index ${SAMPLE}_sorted.bam

mkdir -p ${OUTDIR}
mv ${SAMPLE}_sorted.bam ${SAMPLE}_sorted.bam.bai ${OUTDIR}/
