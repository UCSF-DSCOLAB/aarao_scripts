#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
FQ1       : A path to FQ1
SAMPLE    : The name of the sample
OUTDIR    : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
FQ2       : A path to FQ2
BWAREF    : A path to a bwa reference (default=/krummellab/data1/ipi/data/refs/bwa/hg38)
RGID      : Read Group ID (default=PARSED_FROM_1ST_READ to be FLOWCELL_ID.LANE_NUM)
RGPL      : Read Group Platform (default=ILLUMINA)
RGSM      : Read Group Sample (default=<SAMPLE>)
RGPU      : Read Group Platform Unit (default=PARSED_FROM_1ST_READ to be FLOWCELL_ID.LANE_NUM)
RGBC      : Read Group Barcode (default=PARSED_FROM_1ST_1000_READS and warns if < 75% are concordant)
EOF

if [ ${FQ2-"EMPTY"} == "EMPTY" ]
then
    FQ2="EMPTY"
fi

if [ ${BWAREF-"EMPTY"} == "EMPTY" ]
then
    BWAREF=/krummellab/data1/ipi/data/refs/bwa/hg38
else
    BWAREF=${BWAREF%.fa}
fi

if [ ${RGID-"EMPTY"} == "EMPTY" ]
then
    RGID="EMPTY"
fi

if [ ${RGPL-"EMPTY"} == "EMPTY" ]
then
    RGPL="EMPTY"
fi

if [ ${RGSM-"EMPTY"} == "EMPTY" ]
then
    RGSM="EMPTY"
fi

if [ ${RGPU-"EMPTY"} == "EMPTY" ]
then
    RGPU="EMPTY"
fi

if [ ${RGBC-"EMPTY"} == "EMPTY" ]
then
    RGBC="EMPTY"
fi



echo "Received the following options:"
echo "FQ1      : "${FQ1-""}
echo "FQ2      : "${FQ2-""}
echo "SAMPLE   : "${SAMPLE-""}
echo "OUTDIR   : "${OUTDIR-""}
echo "RGID     : "${RGID-""}
echo "RGPL     : "${RGPL-""}
echo "RGSM     : "${RGSM-""}
echo "RGPU     : "${RGPU-""}
echo "RGBC     : "${RGBC-""}

if [ ${FQ1-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "BWAREF   : "${BWAREF}
echo "LOGDIR   : "${LOGDIR}
echo "NODEREQS : "${NODEREQS}
echo "MEMREQS  : "${MEMREQS}
echo -e "\n"

qsub -v "FQ1=$(readlink -e ${FQ1}),\
         FQ2=$(readlink -f ${FQ2}),\
         SAMPLE=${SAMPLE},\
         OUTDIR=$(readlink -f ${OUTDIR}),\
         BWAREF=$(readlink -f ${BWAREF}),\
         RGID=${RGID},\
         RGPL=${RGPL},\
         RGSM=${RGSM},\
         RGPU=${RGPU},\
         RGBC=${RGBC}" \
     -e ${LOGDIR}/demuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/demuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N bwa_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh

