#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE   : A path to the aligned sam/bam/cram file
SAMPLE    : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
GATKVERSION   : The version of gatk to use (4.0.2.1)
GENOMEREF      : A path to the reference genome fasta (/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa)
EOF

if [ ${GATKVERSION-"EMPTY"} == "EMPTY" ]
then
    GATKVERSION="4.0.2.1"
fi

if [ ${GENOMEREF-"EMPTY"} == "EMPTY" ]
then
    GENOMEREF="/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa"
fi


echo "Received the following options:"
echo "SAMPLE   : "${SAMPLE-""}
echo "SAMFILE  : "${SAMFILE-""}
echo "GATKVERSION  : "${GATKVERSION-""}
echo "GENOMEREF  : "${GENOMEREF-""}

if [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${SAMFILE-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR   : "${LOGDIR}
echo "NODEREQS : "${NODEREQS}
echo "MEMREQS  : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

qsub -v "SAMFILE=$(readlink -f ${SAMFILE}),SAMPLE=${SAMPLE},GENOMEREF=$(readlink -f ${GENOMEREF}),MEMORY=${MEMORY},GATKVERSION=${GATKVERSION}" \
     -e ${LOGDIR}/IndelRealignment_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/IndelRealignment_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N IndelRealignment_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
