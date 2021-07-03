#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh
source $(dirname ${0})/../../../bash/essentials.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
INPUT_DIR          : Path to a directory containing files to process.
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
NODEREQS           : Default node requirements (default: nodes=1:ppn=1)[overrides global default]
MEMREQS            : Default mem requirements (default: vmem=20gb)[overrides global default]
EOF

if [ ${INPUT_DIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if ! [ -d ${INPUT_DIR} ]
then
    echo -e "\nERROR: INPUT_DIR must be a directory\n"
    print_help
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=20gb"
fi

echo "Received the following options:"
echo "INPUT_DIR     : "${INPUT_DIR-""}

echo "LOGDIR         : "${LOGDIR}
echo "NODEREQS       : "${NODEREQS}
echo "MEMREQS        : "${MEMREQS}
echo -e "\n"

export_vars="\
INPUT_DIR=$(readlink -e ${INPUT_DIR})"

UUID=`randomstr 10`

qsub -v ${export_vars} \
     -e ${LOGDIR}/md5sum_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/md5sum_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N md5sum_${UUID} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
