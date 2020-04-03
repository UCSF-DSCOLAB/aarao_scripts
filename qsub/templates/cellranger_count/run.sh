#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_count_${SAMPLE} && cd /scratch/arrao/cellranger_count_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/cellranger_count_${SAMPLE} ; }" EXIT

featureref_argstring=" "
if grep -q "Antibody Capture" ${LIBRARIES_CSV}
then
    featureref_argstring="--feature-ref ${FEATURE_REF} "
fi

cellranger count --id=${SAMPLE} \
                 ${featureref_argstring} \
                 --libraries ${LIBRARIES_CSV} \
                 --chemistry=${CHEMISTRY} \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${MEMORY}

mv /scratch/arrao/cellranger_count_${SAMPLE}/${SAMPLE}/outs/* ${OUTDIR}/

