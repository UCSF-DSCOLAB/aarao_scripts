#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}


mkdir /scratch/arrao/VariantFiltration_${SAMPLE} && cd /scratch/arrao/VariantFiltration_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/VariantFiltration_${SAMPLE} /scratch/arrao/${SAMPLE}_javatmp ; }" EXIT

if [ ${INTERVALFILE-"EMPTY"} == "NONE" ]
then
    INTERVALSTRING=" "
else
    INTERVALFILE2=$(readlink -f ${INTERVALFILE})
    INTERVALSTRING="--intervals ${INTERVALFILE2}"
fi

extension=${VCFFILE##*.}

out_dir=$(dirname ${VCFFILE})
out_base=$(basename ${VCFFILE%.${extension}})_filtered

gatk --java-options "-Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
    VariantFiltration \
    ${INTERVALSTRING} \
    -R ${GENOMEREF} \
    -V {VCFFILE} \
    -O ${out_base}.vcf.gz \
    --cluster-window-size 30 \
    -cluster 3 \
    --filter-name pval \
    -filter "FS > ${PHREDSCOREDPVAL}" 

# TODO: Think about other filters to add

mv /scratch/arrao/VariantFiltration_${SAMPLE}/${out_base}* ${out_dir}/
