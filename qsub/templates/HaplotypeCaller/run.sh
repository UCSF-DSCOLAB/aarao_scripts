#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}


mkdir /scratch/arrao/HaplotypeCaller_${SAMPLE} && cd /scratch/arrao/HaplotypeCaller_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/HaplotypeCaller_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT

out_dir=$(dirname ${SAMFILE})

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   HaplotypeCaller \
   -R ${GENOMEREF} \
   -I ${SAMFILE} \
   -O ${SAMPLE}_GL.vcf.gz \
   -ERC GVCF

mv /scratch/arrao/HaplotypeCaller_${SAMPLE}/${SAMPLE}_GL* ${out_dir}/
