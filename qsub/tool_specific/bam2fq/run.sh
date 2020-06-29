#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/arrao/bam2fq_${SAMPLE} && cd /scratch/arrao/bam2fq_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/bam2fq_${SAMPLE} ; }" EXIT

cp ${BAMFILE} unmapped.bam
if [ -f ${BAMFILE}.bai ]
then
  cp ${BAMFILE}.bai unmapped.bam.bai
elif [ -f ${BAMFILE%.*}.bai ]
then
  cp ${BAMFILE%.*}.bai unmapped.bam.bai
else
  echo "Could not find the bam index for ${BAMFILE} (Tried \"${BAMFILE}.bai\" and \"${BAMFILE%.*}.bai\")." >&2
  exit 1
fi

java -XX:ParallelGCThreads=${PBS_NUM_PPN} \
    -Xmx100G -jar /data/shared/krummellab/ipi/software/picard-2.18.14/picard.jar SamToFastq I=unmapped.bam \
    F=unmapped_1.fq \
    F2=unmapped_2.fq \
    FU=unmapped_up.fq \
    VALIDATION_STRINGENCY=LENIENT
gzip unmapped_1.fq
gzip unmapped_2.fq
gzip unmapped_up.fq

mkdir -p ${OUTDIR}/${SAMPLE} && mv *.fq.gz ${OUTDIR}/${SAMPLE}/
