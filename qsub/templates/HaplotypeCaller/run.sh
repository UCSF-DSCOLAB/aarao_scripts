#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}


mkdir /scratch/arrao/HaplotypeCaller_${SAMPLE} && cd /scratch/arrao/HaplotypeCaller_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/HaplotypeCaller_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT



if [ ${INTERVALFILE-"EMPTY"} == "NONE" ]
then
    INTERVALSTRING=" "
else
    INTERVALFILE2=$(readlink -f ${INTERVALFILE})
    INTERVALSTRING="--intervals ${INTERVALFILE2}"
fi

if [ ${DBSNPFILE-"EMPTY"} == "NONE" ]
then
    DBSNPSTRING=" "
else
    DBSNPFILE2=$(readlink -f ${DBSNPFILE})
    DBSNPSTRING="--dbsnp ${DBSNPFILE2}"
fi

out_dir=$(dirname ${SAMFILE})
gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   HaplotypeCaller \
   ${INTERVALSTRING} \
   ${DBSNPSTRING} \
   -R ${GENOMEREF} \
   -I ${SAMFILE} \
   -O ${SAMPLE}_GL.vcf.gz \
   -ERC GVCF \
   --gvcf-gq-bands 10 \
   --gvcf-gq-bands 20 \
   --gvcf-gq-bands 30 \
   --gvcf-gq-bands 40 \
   --gvcf-gq-bands 60 \
   --gvcf-gq-bands 80 

mv /scratch/arrao/HaplotypeCaller_${SAMPLE}/${SAMPLE}_GL* ${out_dir}/
