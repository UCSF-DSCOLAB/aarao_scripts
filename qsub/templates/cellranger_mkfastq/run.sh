#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_mkfastq_${UUID} && cd /scratch/arrao/cellranger_mkfastq_${UUID} 
trap "{ rm -rf /scratch/arrao/cellranger_mkfastq_${UUID} ; }" EXIT

run_mem=`echo "${MEMORY} * 0.95 / 1" | bc`

cellranger mkfastq --csv=${SAMPLESHEET} \
                   --localcores=${PBS_NUM_PPN} \
                   --run=${BCLDIR} \
                   --lanes=1,2 \
                   --output-dir=/scratch/arrao/cellranger_mkfastq_${UUID} \
                   --localcores=${PBS_NUM_PPN} \
                   --localmem=${run_mem}

mv /scratch/arrao/cellranger_mkfastq_${UUID}/ ${OUTDIR}/
