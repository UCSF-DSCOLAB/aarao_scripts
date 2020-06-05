#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE   : A path to the aligned sam/bam/cram file
SAMPLE    : The name of the sample
VCF       : A path to the vcf file with genotypes for each sample in the aligned sam/bam/cram
OUTDIR    : A folder within which we will place the output dir
EOF

BWAREF=/data/shared/krummellab/ipi/data/refs/bwa/hg38
RG='@RG\tID:foo\tSM:bar'

echo "Received the following options:"
echo "SAMPLE   : "${SAMPLE-""}
echo "SAMFILE  : "${SAMFILE-""}
echo "VCF      : "${VCF-""}
echo "OUTDIR   : "${OUTDIR-""}

if [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${SAMFILE-"ERR"} == "ERR" ] || [ ${VCF-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR   : "${LOGDIR}
echo "NODEREQS : "${NODEREQS}
echo "MEMREQS  : "${MEMREQS}
echo -e "\n"

exit 1
qsub -v "SAMFILE=$(readlink -f ${SAMFILE}),OUTDIR=$(readlink -f ${OUTDIR}),SAMPLE=${SAMPLE},VCF=$(readlink -f ${VCF})" \
     -e ${LOGDIR}/demuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/demuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N demuxlet_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
