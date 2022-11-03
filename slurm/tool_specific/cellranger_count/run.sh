#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

mkdir ${TMPDIR}/cellranger_count_${SAMPLE}
cd ${TMPDIR}/cellranger_count_${SAMPLE}

featureref_argstring=" "
featureref_dirname=" "
if grep -q "Antibody Capture" ${LIBRARIES_CSV}
then
    featureref_argstring="--feature-ref ${FEATUREREF} "
    featureref_dirname=$(dirname ${FEATUREREF})
fi

csv_folders=($(tail -n+2 ${LIBRARIES_CSV} | cut -f 1 -d "," | xargs))

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${LIBRARIES_CSV}) ${csv_folders[@]} ${featureref_dirname} ${TRANSCRIPTOME} ${PWD})

if [[ INCLUDE_INTRONS == 'FALSE' ]]
then
    INTRONS_STRING=""
else
    INTRONS_STRING="--include-introns"
fi



singularity exec \
            ${bindmount_string} \
            -B ${TMPDIR}:/tmp/ \
            --pwd ${PWD} \
            ${CONTAINER} cellranger count \
                                --id=${SAMPLE} \
                                ${featureref_argstring} \
                                --libraries=${LIBRARIES_CSV} \
                                --chemistry=${CHEMISTRY} \
                                --transcriptome=${TRANSCRIPTOME} \
                                --localcores=${SLURM_CPUS_PER_TASK} \
                                --localmem=${USABLEMEMORY} ${INTRONS_STRING}

if [[ -d ${OUTDIR} ]]
then
    mv ${TMPDIR}/cellranger_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_counts
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${TMPDIR}/cellranger_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
fi
