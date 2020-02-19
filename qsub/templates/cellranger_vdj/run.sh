#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_vdj_${SAMPLEID}_${UUID} && cd /scratch/arrao/cellranger_vdj_${SAMPLEID}_${UUID} 
trap "{ rm -rf /scratch/arrao/cellranger_vdj_${SAMPLEID}_${UUID} ; }" EXIT

run_mem=`echo "${MEMORY} * 0.95 / 1" | bc`

cellranger vdj --id=${SAMPLEID} \
               --fastqs=${FASTQDIR} \
               --reference=${REFERENCE} \
               --chain=${CHAIN} \
               --localcores=${PBS_NUM_PPN} \
               --localmem=${run_mem}

mv ${SAMPLEID}/ ${OUTDIR}/
