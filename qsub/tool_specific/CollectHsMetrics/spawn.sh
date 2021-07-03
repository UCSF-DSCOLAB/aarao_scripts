#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE       : A path to the aligned sam/bam/cram file
SAMPLE        : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
BAIT_INTERVALS    : An intervals file to calculate over (/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list)
TARGET_INTERVALS  : An intervals file to calculate over (/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list)
PICARDVERSION       : The version of gatk to use (4.0.2.1)
EOF

if [ ${BAIT_INTERVALS-"EMPTY"} == "EMPTY" ]
then
    BAIT_INTERVALS="/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list"
fi

if [ ${TARGET_INTERVALS-"EMPTY"} == "EMPTY" ]
then
    TARGET_INTERVALS="/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list"
fi

if [ ${PICARDVERSION-"EMPTY"} == "EMPTY" ]
then
    PICARDVERSION="2.22.0"
fi

echo "Received the following options:"
echo "SAMPLE        : "${SAMPLE-""}
echo "SAMFILE       : "${SAMFILE-""}
echo "BAIT_INTERVALS  : "${BAIT_INTERVALS-""}
echo "TARGET_INTERVALS  : "${TARGET_INTERVALS-""}
echo "PICARDVERSION   : "${PICARDVERSION-""}

if [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${SAMFILE-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR        : "${LOGDIR}
echo "NODEREQS      : "${NODEREQS}
echo "MEMREQS       : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
SAMFILE=$(readlink -f ${SAMFILE}),\
SAMPLE=${SAMPLE},\
BAIT_INTERVALS=${BAIT_INTERVALS},\
TARGET_INTERVALS=${TARGET_INTERVALS},\
MEMORY=${MEMORY},\
PICARDVERSION=${PICARDVERSION}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/CollectHsMetrics_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/CollectHsMetrics_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N CollectHsMetrics_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
