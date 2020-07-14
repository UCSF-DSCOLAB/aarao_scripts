#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/${USER}/aarao_scripts/bash/essentials.sh
uuid=`randomstr 10`

source /krummellab/data1/ipi/software/cellranger/usr/SOURCE_THIS

mkdir /scratch/${USER}/cellranger_vdj_${SAMPLE}_${uuid}
cd /scratch/${USER}/cellranger_vdj_${SAMPLE}_${uuid} 
trap "{ rm -rf /scratch/${USER}/cellranger_vdj_${SAMPLE}_${uuid} ; }" EXIT

echo "running command: "
echo "cellranger-${CELLRANGERVERSION} vdj --id=${SAMPLE} \
               --fastqs=${FASTQDIR} \
               --reference=${REFERENCE} \
               --chain=${CHAIN} \
               --localcores=${PBS_NUM_PPN} \
               --localmem=${MEMORY}"
               
cellranger-${CELLRANGERVERSION} vdj --id=${SAMPLE} \
               --fastqs=${FASTQDIR} \
               --reference=${REFERENCE} \
               --chain=${CHAIN} \
               --localcores=${PBS_NUM_PPN} \
               --localmem=${MEMORY}

if [ -f ${OUTDIR} ]
then
    mv /scratch/${USER}/cellranger_vdj_${SAMPLE}_${uuid}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_vdj_${uuid}
else
    mkdir -p $(dirname ${OUTDIR})
    mv /scratch/${USER}/cellranger_vdj_${SAMPLE}_${uuid}/${SAMPLE}/outs ${OUTDIR}
fi
