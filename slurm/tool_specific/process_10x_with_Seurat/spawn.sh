#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh
source $(dirname ${0})/../../../bash/essentials.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMPLE_YML         : Path to a list of samples (one per line) expected to be in DATADIR
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
CONTAINER          : The singularity container to use (/krummellab/data1/singularity_images/RSingleCell/v2/RSingleCell.sif)
WORKING_FOLDER     : A folder in which we will place the run warnings(\${TMPDIR})
RSCRIPTS_DIR       : Where to find aux R scripts form the aarao_scripts repo (\${SCRIPTS}/R)
GENESET_DIR        : Where to find genesets (/krummellab/data1/ipi/data/refs/10x/genesets)
EOF

read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
TIME               : Max Runtime for job (default=05:00:00)
CPUSPERTASK        : Default node requirements (default: 1)
MEMPERCPU          : Default mem requirements (default: 75gb)
EOF

if [ ${SAMPLE_YML-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help
fi

# Override defaults only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TIME" ]]
then
    TIME="05:00:00"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CPUSPERTASK" ]]
then
    CPUSPERTASK="1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMPERCPU" ]]
then
    MEMPERCPU="75gb"
fi

#process optional args
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CONTAINER" ]]
then
    CONTAINER="/krummellab/data1/singularity_images/RStudioSingleCell/v2/RStudioSingleCell.sif"
else
    if [ ! -f ${CONTAINER} ]
    then
       echo -e "\nERROR: CONTAINER does not exist: ${CONTAINER}\n"
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "WORKING_FOLDER" ]]
then
    WORKING_FOLDER="TMPDIR"
else
    if [ ! -d ${WORKING_FOLDER} ]
    then
       echo -e "\nERROR: WORKING_FOLDER does not exist: ${WORKING_FOLDER}\n"
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "RSCRIPTS_DIR" ]]
then
    RSCRIPTS_DIR="${SCRIPTS}/R"
else
    if [ ! -d ${RSCRIPTS_DIR} ]
    then
       echo -e "\nERROR: RSCRIPTS_DIR does not exist: ${RSCRIPTS_DIR}\n"
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GENESET_DIR" ]]
then
    GENESET_DIR="/krummellab/data1/ipi/data/refs/10x/genesets"
else
    if [ ! -d ${GENESET_DIR} ]
    then
       echo -e "\nERROR: GENESET_DIR does not exist: ${GENESET_DIR}\n"
    fi
fi

LOCAL_EXPORT_VARS="\
CONTAINER=$(readlink -e ${CONTAINER}),\
SAMPLE_YML=$(readlink -e ${SAMPLE_YML}),\
WORKING_FOLDER=${WORKING_FOLDER},\
RSCRIPTS_DIR=$(readlink -e ${RSCRIPTS_DIR}),\
GENESET_DIR=$(readlink -e ${GENESET_DIR})"

JOBNAME=process_10x_seurat_$(randomstr 10)

source $(dirname ${0})/../final_spawn.sh
