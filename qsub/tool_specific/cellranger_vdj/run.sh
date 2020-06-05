#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/arrao/scripts/bash/essentials.sh
uuid=`randomstr 10`

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_vdj_${SAMPLE}_${uuid} && cd /scratch/arrao/cellranger_vdj_${SAMPLE}_${uuid} 
trap "{ rm -rf /scratch/arrao/cellranger_vdj_${SAMPLE}_${uuid} ; }" EXIT


cellranger vdj --id=${SAMPLE} \
               --fastqs=${FASTQDIR} \
               --reference=${REFERENCE} \
               --chain=${CHAIN} \
               --localcores=${PBS_NUM_PPN} \
               --localmem=${MEMORY}

if [ -f ${OUTDIR} ]
then
    mv /scratch/arrao/cellranger_vdj_${SAMPLE}_${uuid}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_vdj_${uuid}
else
    mkdir -p $(dirname ${OUTDIR})
    mv /scratch/arrao/cellranger_vdj_${SAMPLE}_${uuid}/${SAMPLE}/outs ${OUTDIR}
fi
