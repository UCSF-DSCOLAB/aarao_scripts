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
WORKING_FOLDER     : A folder in which we will place the run warnings(\${TMPDIR})
RSCRIPTS_DIR       : Where to find aux R scripts form the aarao_scripts repo (\${SCRIPTS}/R)
GENESET_DIR        : Where to find genesets (/krummellab/data1/ipi/data/refs/10x/genesets)
NODEREQS           : Default node requirements (default: nodes=1:ppn=1)[overrides global default]
MEMREQS            : Default mem requirements (default: vmem=150gb)[overrides global default]
EOF

if [ ${SAMPLE_YML-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=150gb"
fi

#process optional args
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

echo "Received the following options:"
echo "SAMPLE_YML     : "${SAMPLE_YML-""}

echo "RSCRIPTS_DIR   : "${RSCRIPTS_DIR-""}
echo "GENESET_DIR    : "${GENESET_DIR-""}

echo "LOGDIR         : "${LOGDIR}
echo "NODEREQS       : "${NODEREQS}
echo "MEMREQS        : "${MEMREQS}
echo -e "\n"

export_vars="\
SAMPLE_YML=$(readlink -e ${SAMPLE_YML}),\
WORKING_FOLDER=$(readlink -e ${WORKING_FOLDER}),\
RSCRIPTS_DIR=$(readlink -e ${RSCRIPTS_DIR}),\
GENESET_DIR=$(readlink -e ${GENESET_DIR})"

UUID=`randomstr 10`

qsub -v ${export_vars} \
     -e ${LOGDIR}/process_10x_seurat_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/process_10x_seurat_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N process_10x_seurat_${UUID} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
