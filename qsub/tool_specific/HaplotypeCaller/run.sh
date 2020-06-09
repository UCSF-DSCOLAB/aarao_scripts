#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}


mkdir /scratch/${USER}/HaplotypeCaller_${SAMPLE} && cd /scratch/${USER}/HaplotypeCaller_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/HaplotypeCaller_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT



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

if [ ${GVCF-"FALSE"} == "FALSE" ]
then
    GVCFSTRING="-ERC NONE"
else
    # Can't handle ERC BP_RESOLUTION right now.
    GVCFSTRING="\
-ERC GVCF \
--gvcf-gq-bands 10 \
--gvcf-gq-bands 20 \
--gvcf-gq-bands 30 \
--gvcf-gq-bands 40 \
--gvcf-gq-bands 60 \
--gvcf-gq-bands 80"

fi
    

out_dir=$(dirname ${SAMFILE})
gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   HaplotypeCaller \
   ${INTERVALSTRING} \
   ${GVCFSTRING} \
   ${DBSNPSTRING} \
   -R ${GENOMEREF} \
   -I ${SAMFILE} \
   -O ${SAMPLE}_GL.vcf.gz \
   -stand-call-conf ${STANDCALLCONF} \
   --dont-use-soft-clipped-bases ${DONTUSESOFTCLIPPEDBASES}
   

mv /scratch/${USER}/HaplotypeCaller_${SAMPLE}/${SAMPLE}_GL* ${out_dir}/
