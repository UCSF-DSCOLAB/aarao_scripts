#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/ipi/software/cellranger/usr/SOURCE_THIS

mkdir /scratch/${USER}/cellranger_count_${SAMPLE}
cd /scratch/${USER}/cellranger_count_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/cellranger_count_${SAMPLE} ; }" EXIT

featureref_argstring=" "
if grep -q "Antibody Capture" ${LIBRARIES_CSV}
then
    featureref_argstring="--feature-ref ${FEATUREREF} "
fi


echo "running command: "
echo "cellranger-${CELLRANGERVERSION} count --id=${SAMPLE} \
                 ${featureref_argstring} \
                 --libraries ${LIBRARIES_CSV} \
                 --chemistry=${CHEMISTRY} \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${MEMORY}"
                 
cellranger-${CELLRANGERVERSION} count --id=${SAMPLE} \
                 ${featureref_argstring} \
                 --libraries ${LIBRARIES_CSV} \
                 --chemistry=${CHEMISTRY} \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${MEMORY}

if [ -f ${OUTDIR} ]
then
    mv /scratch/${USER}/cellranger_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_counts
else
    mkdir -p $(dirname ${OUTDIR})
    mv /scratch/${USER}/cellranger_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
fi
