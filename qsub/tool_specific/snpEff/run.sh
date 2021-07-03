#!/bin/bash
set -e
set -o nounset

module load CBC jdk/8

mkdir /scratch/${USER}/snpEff_${PATIENT}
cd /scratch/${USER}/snpEff_${PATIENT} 
trap "{ rm -rf /scratch/${USER}/snpEff_${PATIENT} ; }" EXIT

out_dir=$(dirname ${INVCF})
out_base=$(basename ${INVCF%.vcf})_snpEffed

cp ${INVCF} .
java -Djava.io.tmpdir=/scratch/${USER}/${PATIENT}_javatmp \
     -Xmx${MEMORY}g \
     -jar /krummellab/data1/ipi/software/snpeff/${SNPEFFVERSION}/snpEff/snpEff.jar \
          eff \
          -dataDir ${SNPEFFDBDIRECTORY} \
          -c ${SNPEFFDBCONFIG} \
          -no-intergenic \
          -no-downstream \
          -no-upstream \
          -noStats \
          ${SNPEFFDB} \
          $(basename ${INVCF}) > ${out_base}.vcf

rm -f $(basename ${INVCF})
mv /scratch/${USER}/snpEff_${PATIENT}/${out_base}* ${out_dir}/

