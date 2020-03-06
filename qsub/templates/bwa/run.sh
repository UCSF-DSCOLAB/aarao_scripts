#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/ipi/ipi_usr/SOURCE_THIS 

mkdir /scratch/arrao/bwa_${SAMPLE} && cd /scratch/arrao/bwa_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/bwa_${SAMPLE} ; }" EXIT

getRG(){
 zcat ${1} | 
 head -1 | 
 awk -F":" -v SAMPLE=${2} '{print "@RG\\tID:"$3"."$2"."SAMPLE"."$4"\\tSM:"SAMPLE"\\tPU:"$3"."$2"."SAMPLE"."$4"\\tBC:"$NF}'
}

RG=`getRG ${FQ1} ${SAMPLE}`

bwa mem -t ${PBS_NUM_PPN} \
        -v 1 \
        -R ${RG} \
        ${BWAREF} \
        ${FQ1} \
        ${FQ2} > ${SAMPLE}.sam

samtools view -C \
              -o ${SAMPLE}.cram \
              -@ ${PBS_NUM_PPN} \
              -T ${BWAREF}.fa \
              ${SAMPLE}.sam \
              && rm ${SAMPLE}.sam

samtools sort ${SAMPLE}.cram \
              -@ ${PBS_NUM_PPN} \
              --reference ${BWAREF}.fa \
              -o ${SAMPLE}_sorted.cram \
              && rm ${SAMPLE}.cram

samtools index ${SAMPLE}_sorted.cram

mkdir -p ${OUTDIR}
mv ${SAMPLE}_sorted.cram ${SAMPLE}_sorted.cram.crai ${OUTDIR}/
