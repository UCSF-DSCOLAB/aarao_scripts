#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
FASTQDIR           : Path to the fastq directory
SAMPLE             : The name of the sample
OUTDIR             : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
CELLRANGERVERSION  : The version of cellranger to use (3.0.2)
NODEREQS           : Default node requirements (default: nodes=1:ppn=32)[overrides global default]
MEMREQS            : Default mem requirements (default: vmem=250gb)[overrides global default]
REFERENCE          : Path to 10X Indexes (default: /krummellab/data1/ipi/data/refs/10x/GRCh38_VDJ)
CHAIN              : Chain to process [TR/IG/auto/all] (default: auto)
EOF

if [ ${FASTQDIR-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CELLRANGERVERSION" ]]
then
    CELLRANGERVERSION="3.0.2"
else
    allowed=("3.0.2" "3.1.0" "4.0.0")
    if [[ ! "${allowed[@]}" =~ "${CELLRANGERVERSION}" ]]
    then
        echo -e "\nERROR: CELLRANGERVERSION must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=32"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=250gb"
fi


if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "REFERENCE" ]]
then
    REFERENCE=/krummellab/data1/ipi/data/refs/10x/GRCh38_VDJ
else
    if [ ! -f ${REFERENCE} ]
    then
       echo -e "\nERROR: REFERENCE file does not exist: ${REFERENCE}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CHAIN" ]]
then
    CHAIN="auto"
else
    allowed=('TR' 'IG' 'all' 'auto')
    if [[ ! "${allowed[@]}" =~ "${CHAIN}" ]]
    then
        echo -e "\nERROR: CHAIN must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
fi 


echo "Received the following options:"
echo "FASTQDIR           : "${FASTQDIR-""}
echo "OUTDIR             : "${OUTDIR-""}
echo "CELLRANGERVERSION  : "${CELLRANGERVERSION-""}
echo "SAMPLE             : "${SAMPLE-""}
echo "REFERENCE          : "${REFERENCE-""}
echo "CHAIN              : "${CHAIN-""}

echo "LOGDIR             : "${LOGDIR}
echo "NODEREQS           : "${NODEREQS}
echo "MEMREQS            : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
FASTQDIR=$(readlink -e ${FASTQDIR}),\
SAMPLE=${SAMPLE},\
OUTDIR=${OUTDIR},\
CELLRANGERVERSION=${CELLRANGERVERSION},\
REFERENCE=$(readlink -e ${REFERENCE}),\
CHAIN=${CHAIN},\
MEMORY=${MEMORY}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/cellranger_vdj_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/cellranger_vdj_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_vdj_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
