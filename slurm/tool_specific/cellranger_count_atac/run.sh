#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

mkdir ${TMPDIR}/cellranger_atac_count_${SAMPLE}
cd ${TMPDIR}/cellranger_atac_count_${SAMPLE}

if [[ $(grep -c "ATAC Seq" ${LIBRARIES_CSV}) -ne 1 ]]
then
  echo "ERROR: ${LIBRARIES_CSV} must contain ONLY one library of type 'ATAC Seq'"
  exit 1
fi

fastq_dir=$(grep "ATAC Seq" ${LIBRARIES_CSV} | cut -f1 -d ",")
bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${LIBRARIES_CSV}) ${fastq_dir} ${REFERNCE} ${PWD})



singularity exec \
            ${bindmount_string} \
            -B ${TMPDIR}:/tmp/ \
            --pwd ${PWD} \
            ${CONTAINER} cellranger-atac-2.0.0 count \
                                                --id=${SAMPLE} \
                                                --fastqs=${fastq_dir} \
                                                --sample=${SAMPLE} \
                                                --reference=${REFERENCE} \
                                                --localcores=32 \
                                                --localmem=150

if [ -f ${OUTDIR} ]
then
    mv ${TMPDIR}/cellranger_atac_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}/${SAMPLE}_counts
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${TMPDIR}/cellranger_atac_count_${SAMPLE}/${SAMPLE}/outs ${OUTDIR}
fi

