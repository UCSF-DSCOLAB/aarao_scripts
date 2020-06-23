#!/bin/bash
source /krummellab/data1/ipi/software/popscle/usr/SOURCE_THIS
source /krummellab/data1/ipi/software/bedtools/2.29.2/SOURCE_THIS

set -e
set -o nounset

mkdir -p /scratch/${USER}/freemuxlet_${SAMPLE} && cd /scratch/${USER}/freemuxlet_${SAMPLE}
trap "{ rm -rf /scratch/${USER}/freemuxlet_${SAMPLE} ; }" EXIT

prefix=$(basename ${BAMFILE%.*})
extension=${BAMFILE##*.}

# Create filtered BAM with only the reads dsc-pileup needs.
echo "running command: "
echo "/krummellab/data1/ipi/software/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    ${BAMFILE} \
    ${BARCODELIST} \
    ${ONEKGENOMESVCF} \
    ${prefix}_filtered.${extension}"

/krummellab/data1/ipi/software/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    ${BAMFILE} \
    ${BARCODELIST} \
    ${ONEKGENOMESVCF} \
    ${prefix}_filtered.${extension}

# Use filtered BAM file for dsc-pileup.
echo "running command: "
echo "popscle dsc-pileup --sam ${prefix}_filtered.${extension} \
                   --tag-group CB \
                   --tag-UMI UB \
                   --vcf ${ONEKGENOMESVCF} \
                   --group-list ${BARCODELIST} \
                   --out ${SAMPLE}"

popscle dsc-pileup --sam ${prefix}_filtered.${extension} \
                   --tag-group CB \
                   --tag-UMI UB \
                   --vcf ${ONEKGENOMESVCF} \
                   --group-list ${BARCODELIST} \
                   --out ${SAMPLE}

popscle freemuxlet --plp ${SAMPLE} \
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
