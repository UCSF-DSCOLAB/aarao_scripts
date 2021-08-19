#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

mkdir ${TMPDIR}/freemuxlet_${SAMPLE}
cd ${TMPDIR}/freemuxlet_${SAMPLE}

prefix=$(basename ${BAMFILE%.*})
extension=${BAMFILE##*.}

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${BAMFILE}) $(dirname ${BARCODELIST}) $(dirname ${ONEKGENOMESVCF}) ${POPSCLE_HELPER_TOOLS_DIR} ${PWD})

# Create filtered BAM with only the reads dsc-pileup needs.
singularity exec \
            ${bindmount_string} \
            --pwd ${PWD} \
            ${CONTAINER} bash ${POPSCLE_HELPER_TOOLS_DIR}/filter_bam_file_for_popscle_dsc_pileup.sh \
                ${BAMFILE} \
                ${BARCODELIST} \
                ${ONEKGENOMESVCF} \
                ${prefix}_filtered.${extension}

if [[ NO_TAG_UMI == 'FALSE' ]]
then
    UMI_TAG_STRING="--tag-UMI UB"
else
    UMI_TAG_STRING=" "
fi

# Use filtered BAM file for dsc-pileup
singularity exec \
            ${bindmount_string} \
            --pwd ${PWD} \
            ${CONTAINER} popscle dsc-pileup \
                    --sam ${prefix}_filtered.${extension} \
                    --tag-group CB \
                    ${UMI_TAG_STRING} \
                    --vcf ${ONEKGENOMESVCF} \
                    --group-list ${BARCODELIST} \
                    --out ${SAMPLE}

singularity exec \
            ${bindmount_string} \
            --pwd ${PWD} \
            ${CONTAINER} popscle freemuxlet \
                    --plp ${SAMPLE} \
                    --out ${SAMPLE} \
                    --nsample ${NUMSAMPLES} \
                    --seed ${RANDOMSEED} \
                    --group-list ${BARCODELIST}

if [ -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}/${SAMPLE}_freemuxlet && \
      mv ${SAMPLE}* ${OUTDIR}/${SAMPLE}_freemuxlet
else
    mkdir -p ${OUTDIR} && \
      mv ${SAMPLE}* ${OUTDIR}/
fi
