#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

mkdir ${TMPDIR}/cellranger_vdj_${SAMPLE}
cd ${TMPDIR}/cellranger_vdj_${SAMPLE}

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${FASTQDIR}) ${REFERENCE} ${PWD})

singularity exec \
            ${bindmount_string} \
            --pwd ${PWD} \
            -B ${TMPDIR}:/tmp/ \
            ${CONTAINER} cellranger vdj \
                                --id=${SAMPLE} \
                                --fastqs=${FASTQDIR} \
                                --reference=${REFERENCE} \
                                --chain=${CHAIN} \
                                --localcores=${SLURM_CPUS_PER_TASK} \
                                --localmem=${USABLEMEMORY}

if [[ -d ${OUTDIR} ]]
then
    mv ${TMPDIR}/cellranger_vdj_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_vdj
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${TMPDIR}/cellranger_vdj_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
fi
