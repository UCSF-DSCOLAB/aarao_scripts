#!/bin/bash
module load CBC gatk/${GATKVERSION}
source /krummellab/data1/ipi/software/samtools/samtools-1.10-usr/SOURCE_THIS

set -e
set -o nounset

mkdir /scratch/arrao/MuTect2_${PATIENT}
cd /scratch/arrao/MuTect2_${PATIENT} 
trap "{ rm -rf /scratch/arrao/MuTect2_${PATIENT} ; }" EXIT

cp ${TUMOR_BAMFILE} .
if [ -f ${TUMOR_BAMFILE}.bai ]
then
  cp ${TUMOR_BAMFILE}.bai .
elif [ -f ${TUMOR_BAMFILE%*.bam}.bai ]
then
  cp ${TUMOR_BAMFILE%*.bam}.bai .
else
  samtools index $(basename ${TUMOR_BAMFILE})
fi
cp ${NORMAL_BAMFILE} .
if [ -f ${NORMAL_BAMFILE}.bai ]
then
  cp ${NORMAL_BAMFILE}.bai .
elif [ -f ${NORMAL_BAMFILE%*.bam}.bai ]
then
  cp ${NORMAL_BAMFILE%*.bam}.bai .
else
  samtools index $(basename ${NORMAL_BAMFILE})
fi

gatk --java-options "-Djava.io.tmpdir=${PATIENT}_javatmp -Xmx${MEMORY}g" \
     Mutect2 \
     --input $(basename ${TUMOR_BAMFILE}) \
     --input $(basename ${NORMAL_BAMFILE}) \
     --output ${PATIENT}_mutect2_mutations.vcf \
     --tumor-sample ${TUMOR_SAMPLE_NAME} \
     --normal-sample ${NORMAL_SAMPLE_NAME} \
     --reference ${GENOMEREF} \
     --intervals ${TARGET_INTERVALS} \
     --native-pair-hmm-threads ${PBS_NUM_PPN} \
     --germline-resource ${GERMLINE_RESOURCE} 

if [ ! -d ${OUTDIR} ]
then
     mkdir -p ${OUTDIR}
fi

mv /scratch/arrao/MuTect2_${PATIENT}/${PATIENT}_mutect2_mutations.vcf* ${OUTDIR}/
