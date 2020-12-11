#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
FQ1        : A path to FQ1
SAMPLE     : The name of the sample
OUTDIR     : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
FQ2        : A path to FQ2
ADAPTERSEQ : Sequence used for adapters [NEXTERA, TRUSEQ, or a sequence of nucleotides] (default=<NEXTERA>)
EOF

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "FQ2" ]]
then
    FQ2="EMPTY"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "ADAPTERSEQ" ]]
then
    ADAPTERSEQ="NEXTERA"
elif [ "${ADAPTERSEQ}" != "TRUSEQ "]
then
    if ! [[ ${ADAPTERSEQ} =~ ^[ACTGN]+$ ]]
    then
        echo -e "\nERROR: ADAPTERSEQ must be one of \`TRUSEQ\`, \`NEXTERA\`, or a string composed of ACTGN only \n"
        print_help 
    fi
fi

echo "Received the following options:"
echo "FQ1          : "${FQ1-""}
echo "FQ2          : "${FQ2-""}
echo "SAMPLE       : "${SAMPLE-""}
echo "OUTDIR       : "${OUTDIR-""}
echo "ADAPTERSEQ   : "${ADAPTERSEQ-""}

if [ ${FQ1-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR       : "${LOGDIR}
echo "NODEREQS     : "${NODEREQS}
echo "MEMREQS      : "${MEMREQS}
echo -e "\n"

export_vars="\
FQ1=$(readlink -e ${FQ1}),\
FQ2=$(readlink -f ${FQ2}),\
SAMPLE=${SAMPLE},\
ADAPTERSEQ=${ADAPTERSEQ},\
OUTDIR=$(readlink -f ${OUTDIR})"

qsub -v ${export_vars} \
     -e ${LOGDIR}/fastp_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/fastp_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N fastp_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh

