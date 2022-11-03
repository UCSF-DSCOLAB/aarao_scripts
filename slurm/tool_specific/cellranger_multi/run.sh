#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

mkdir ${TMPDIR}/cellranger_multi_${SAMPLE}
cd ${TMPDIR}/cellranger_multi_${SAMPLE}

# https://unix.stackexchange.com/a/264977/411883
fastq_folders=($(sed -n '/\[libraries\]/,/^$/{p;/^$/q}' ${MULTI_CSV} | 
                    tail -n+3 |
                    grep -v "^$" |
                    cut -f 2 -d "," | xargs))
reference_folders=($(grep "^reference" ${MULTI_CSV} | 
                    cut -f 2 -d "," | xargs))

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${MULTI_CSV}) ${fastq_folders[@]} ${reference_folders[@]} ${PWD})

singularity exec \
            ${bindmount_string} \
            -B ${TMPDIR}:/tmp/ \
            --pwd ${PWD} \
            ${CONTAINER} cellranger multi \
                                --id=${SAMPLE} \
                                --csv=${MULTI_CSV} \
                                --localcores=${SLURM_CPUS_PER_TASK} \
                                --localmem=${USABLEMEMORY}

if [[ -d ${OUTDIR} ]]
then
    mv ${TMPDIR}/cellranger_multi_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_multis
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${TMPDIR}/cellranger_multi_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
fi
