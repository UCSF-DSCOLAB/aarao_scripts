#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_count_${SAMPLEID} && cd /scratch/arrao/cellranger_count_${SAMPLEID} 
trap "{ rm -rf /scratch/arrao/cellranger_count_${SAMPLEID} ; }" EXIT

run_mem=`echo "${MEMORY} * 0.95 / 1" | bc`

cellranger count --id=${SAMPLEID} \
                 --feature-ref ${FEATURE_REF} \
                 --libraries ${LIBRARIES_CSV} \
                 --chemistry=${CHEMISTRY} \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${run_mem}

mv /scratch/arrao/cellranger_count_${SAMPLEID}/${SAMPLEID}/outs/* ${OUTDIR}/

