#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
FQ1       : A path to FQ1
FQ2       : A path to FQ2
SAMPLE    : The name of the sample
OUTDIR    : A folder within which we will place the output dir
EOF

echo "Received the following options:"
echo "FQ1      : "${FQ1-""}
echo "FQ2      : "${FQ2-""}
echo "SAMPLE   : "${SAMPLE-""}
echo "OUTDIR   : "${OUTDIR-""}

if [ ${FQ1-"ERR"} == "ERR" ] || [ ${FQ2-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR   : "${LOGDIR}
echo "NODEREQS : "${NODEREQS}
echo "MEMREQS  : "${MEMREQS}
echo -e "\n"

qsub -v "FQ1=$(readlink -f ${FQ1}),FQ2=$(readlink -f ${FQ2}),SAMPLE=${SAMPLE},OUTDIR=$(readlink -f ${OUTDIR})" \
     -e ${LOGDIR}/fastp_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/fastp_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N fastp_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh

