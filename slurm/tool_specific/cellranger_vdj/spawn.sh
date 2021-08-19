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
CONTAINER          : The singularity container to use (/krummellab/data1/singularity_images/cellranger/6.0.2/cellranger.sif)
REFERENCE          : Path to 10X Indexes (default: /krummellab/data1/ipi/data/refs/10x/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0)
CHAIN              : Chain to process [TR/IG/auto/all] (default: auto)
EOF

read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
TIME               : Max Runtime for job (default=2-00:00:00)
CPUSPERTASK        : Default node requirements (default: 12)
MEMPERCPU          : Default mem requirements (default: 10gb)
EOF


if [ ${FASTQDIR-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if readlink -f ${OUTDIR}
then
    if [[ -d ${OUTDIR} ]] && [[! -w ${OUTDIR} ]]
    then
        echo -e "\nERROR: OUTDIR exists but is not writable."
        exit 1
    elif [[ -d $(dirname ${OUTDIR}) ]] && [[! -w $(dirname ${OUTDIR}) ]]
    then
        echo -e "\nERROR: Parent dir to OUTDIR exists but is not writable."
        exit 1
    fi
else
    echo -e "\nERROR: parent directory for OUTDIR does not exist\n"
    exit 1
fi

# Override defaults only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TIME" ]]
then
    TIME="2-00:00:00"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CPUSPERTASK" ]]
then
    CPUSPERTASK="12"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMPERCPU" ]]
then
    MEMPERCPU="13gb"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "REFERENCE" ]]
then
    REFERENCE=/krummellab/data1/ipi/data/refs/10x/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
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

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CONTAINER" ]]
then
    CONTAINER="/krummellab/data1/singularity_images/cellranger/6.0.2/cellranger.sif"
else
    if [ ! -f ${CONTAINER} ]
    then
       echo -e "\nERROR: CONTAINER does not exist: ${CONTAINER}\n"
    fi
fi

LOCAL_EXPORT_VARS="\
FASTQDIR=$(readlink -e ${FASTQDIR}),\
OUTDIR=$(readlink -f ${OUTDIR}),\
SAMPLE=${SAMPLE},\
CONTAINER=$(readlink -e ${CONTAINER}),\
REFERENCE=$(readlink -e ${REFERENCE}),\
CHAIN=${CHAIN}"

JOBNAME=cellranger_count_${SAMPLE}

source $(dirname ${0})/../final_spawn.sh
