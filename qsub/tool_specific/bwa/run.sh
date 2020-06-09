#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/ipi/ipi_usr/SOURCE_THIS 

mkdir /scratch/${USER}/bwa_${SAMPLE} && cd /scratch/${USER}/bwa_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/bwa_${SAMPLE} ; }" EXIT

# This is old and will eb deleted soon.
getRG(){
 zcat ${1} | 
 awk -F":" -v SAMPLE=${2} '{print "@RG\\tID:"$3"."$2"."SAMPLE"."$4"\\tPL:ILLUMINA\\tSM:"SAMPLE"\\tPU:"$3"."$4"\\tBC:"$NF; exit}'
}

# This will return the 5 RG items we are interested in, in the fixed order
# ID PL SM PU BC
get_RG_items(){
  first_record=$(head -n 1 < <(zcat ${FQ1} 2>/dev/null) | grep "@" )

  BC=($(head -n 4000 < <(zcat ${FQ1} 2>/dev/null) | grep "@" | awk '{print $2}' | sort | uniq -c | sort -k1n | tail -1))

  if [[ ${BC[0]} < 750 ]]
  then
    echo "WARNING: Less than 75% of the top 1000 reads have the most popular barcode (${BC[1]})" >&2
  fi

  echo ${first_record} | 
           awk -F":" -v SAMPLE=${2} -v BARCODE=${BC[1]/*:/} '{print $3"."$4" ILLUMINA "SAMPLE" "$3"."$4" "BARCODE}'
}

RG_items=($(echo $(get_RG_items ${FQ1} ${SAMPLE})))

if [ ${RGID} == "EMPTY" ]
then
    RGID=${RG_items[0]}
fi

if [ ${RGPL} == "EMPTY" ]
then
    RGPL=${RG_items[1]}
fi

if [ ${RGSM} == "EMPTY" ]
then
    RGSM=${RG_items[2]}
fi

if [ ${RGPU} == "EMPTY" ]
then
    RGPU=${RG_items[3]}
fi

if [ ${RGBC} == "EMPTY" ]
then
    RGBC=${RG_items[4]}
fi

RG="@RG\\tID:"$RGID"\\tPL:${RGPL}\\tSM:${RGSM}\\tPU:${RGPU}\\tBC:${RGBC}"


if [ $(basename ${FQ2}) == "EMPTY" ]
then
      FQ2=" "
fi

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
