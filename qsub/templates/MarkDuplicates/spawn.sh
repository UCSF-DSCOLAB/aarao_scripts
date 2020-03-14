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
PICARDVERSION   : The version of picard to use (2.22.0)
GENOMEREF      : A path to the reference genome fasta (/krummellab/data1/ipi/data/refs/bwa/hg38.fa)
EOF

if [ ${PICARDVERSION-"EMPTY"} == "EMPTY" ]
then
    PICARDVERSION="2.22.0"
fi

if [ ${GENOMEREF-"EMPTY"} == "EMPTY" ]
then
    GENOMEREF="/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa"
fi


echo "Received the following options:"
echo "SAMPLE   : "${SAMPLE-""}
echo "SAMFILE  : "${SAMFILE-""}
echo "PICARDVERSION  : "${PICARDVERSION-""}
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

qsub -v "SAMFILE=$(readlink -f ${SAMFILE}),SAMPLE=${SAMPLE},GENOMEREF=$(readlink -f ${GENOMEREF}),MEMORY=${MEMORY},PICARDVERSION=${PICARDVERSION}" \
     -e ${LOGDIR}/MarkDuplicates_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/MarkDuplicates_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N MarkDuplicates_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
