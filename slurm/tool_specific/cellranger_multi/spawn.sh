#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
MULTI_CSV          : Path to a multi config libraries.csv file.
OUTDIR             : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
SAMPLE             : The name of the sample (parsed from MULTI_CSV by default)
EOF

read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
TIME               : Max Runtime for job (default=2-00:00:00)
CPUSPERTASK        : Default node requirements (default: 12)
MEMPERCPU          : Default mem requirements (default: 13gb)
EOF


if [ ${MULTI_CSV-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if readlink -f ${OUTDIR}
then
    if [[ -d ${OUTDIR} ]] && ! [[ -w ${OUTDIR} ]]
    then
        echo -e "\nERROR: OUTDIR exists but is not writable."
        exit 1
    elif [[ -d $(dirname ${OUTDIR}) ]] && ! [[ -w $(dirname ${OUTDIR}) ]]
    then
        echo -e "\nERROR: Parent dir to OUTDIR exists but is not writable."
        exit 1
    fi
else
    echo -e "\nERROR: parent directory for OUTDIR does not exist\n"
    exit 1
fi


if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SAMPLE" ]]
then
    if grep -q "Gene Expression" ${MULTI_CSV}
    then
        SAMPLE=$(grep "Gene Expression" ${MULTI_CSV} | cut -d "," -f 1)
    else
        echo "ERROR: Could not automatically identify SAMPLE, plese specify it explicitly."
        exit 1
    fi
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
MULTI_CSV=$(readlink -e ${MULTI_CSV}),\
OUTDIR=$(readlink -f ${OUTDIR}),\
SAMPLE=${SAMPLE},\
CONTAINER=$(readlink -e ${CONTAINER})"

JOBNAME=cellranger_count_${SAMPLE}

source $(dirname ${0})/../final_spawn.sh
