#!/bin/bash
set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
    then
        TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
	    trap "{ rm -rf ${TMPDIR} ; }" EXIT
	    fi
 
mkdir ${TMPDIR}/demuxlet_${SAMPLE}
cd ${TMPDIR}/demuxlet_${SAMPLE}

echo "BAMFILE: $BAMFILE"
prefix=$(basename ${BAMFILE%.*})
extension=${BAMFILE##*.}

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${BAMFILE}) $(dirname ${BARCODELIST}) $(dirname ${ONEKGENOMESVCF}) $(dirname ${VCF}) ${POPSCLE_HELPER_TOOLS_DIR} ${PWD})
echo $bindmount_string


# Create filtered BAM with only the reads dsc-pileup needs.
if [[ DMX_ONLY == 'FALSE' ]]
then

    singularity exec \
        ${bindmount_string} \
        -B ${TMPDIR}:/tmp/ \
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
        -B ${TMPDIR}:/tmp/ \
        --pwd ${PWD} \
        ${CONTAINER} popscle dsc-pileup \
        --sam ${prefix}_filtered.${extension} \
        --tag-group CB \
        ${UMI_TAG_STRING} \
        --vcf ${ONEKGENOMESVCF} \
        --group-list ${BARCODELIST} \
        --out ${SAMPLE}
fi

if [[ DSC_ONLY == 'FALSE' ]]
then
singularity exec \
                ${bindmount_string} \
                -B ${TMPDIR}:/tmp/ \
                --pwd ${PWD} \
                ${CONTAINER} popscle demuxlet \
                        --plp ${SAMPLE} \
                        --out ${SAMPLE} \
                        --vcf ${VCF} \
                        --field GT \
                        --group-list ${BARCODELIST}
fi

if [ -d ${OUTDIR} ]
    then
        mkdir -p ${OUTDIR}/${SAMPLE}_demuxlet && \
	          mv ${SAMPLE}* ${OUTDIR}/${SAMPLE}_demuxlet
	else
        mkdir -p ${OUTDIR} && \
	          mv ${SAMPLE}* ${OUTDIR}/
	fi
