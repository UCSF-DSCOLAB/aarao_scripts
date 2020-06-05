#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh
source $(dirname ${0})/../../../bash/essentials.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMPLE_YML         : Path to a list of samples (one per line) expected to be in DATADIR
OUT_FOLDER         : A folder within which we will place the output dir(s)
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
RSCRIPTS_DIR       : Where to find aux scripts (/krummellab/data1/arrao/scripts/R)
GENESET_DIR        : Where to find genesets (/krummellab/data1/ipi/data/refs/10x/genesets)
SPECIES            : human or mouse? (human)
RERUN_STAGE        : Should we rerun the last run stage? (FALSE)  
AB_ASSAY_NAME      : What should the antibody capture assay (if any) be called? (IDX)
NODEREQS           : Default node requirements (default: nodes=1:ppn=1)[overrides global default]
MEMREQS            : Default mem requirements (default: vmem=150gb)[overrides global default]
EOF

if [ ${SAMPLE_YML-"ERR"} == "ERR" ] || [ ${OUT_FOLDER-"ERR"} == "ERR" ]
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
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "RSCRIPTS_DIR" ]]
then
    RSCRIPTS_DIR="/krummellab/data1/arrao/scripts/R"
else
    if [ ! -f ${RSCRIPTS_DIR} ]
    then
       echo -e "\nERROR: RSCRIPTS_DIR does not exist: ${RSCRIPTS_DIR}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GENESET_DIR" ]]
then
    GENESET_DIR="/krummellab/data1/ipi/data/refs/10x/genesets"
else
    if [ ! -f ${GENESET_DIR} ]
    then
       echo -e "\nERROR: GENESET_DIR does not exist: ${GENESET_DIR}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SPECIES" ]]
then
    SPECIES="human"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "RERUN_STAGE" ]]
then
    RERUN_STAGE="FALSE"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "AB_ASSAY_NAME" ]]
then
    AB_ASSAY_NAME="IDX"
fi


echo "Received the following options:"
echo "OUT_FOLDER     : "${OUT_FOLDER-""}
echo "SAMPLE_YML        : "${SAMPLE_YML-""}

echo "RSCRIPTS_DIR   : "${RSCRIPTS_DIR-""}
echo "GENESET_DIR    : "${GENESET_DIR-""}
echo "SPECIES        : "${SPECIES-""}
echo "RERUN_STAGE    : "${RERUN_STAGE-""}
echo "AB_ASSAY_NAME  : "${AB_ASSAY_NAME-""}
echo "NODEREQS       : "${NODEREQS-""}
echo "MEMREQS        : "${MEMREQS-""}

echo "LOGDIR         : "${LOGDIR}
echo "NODEREQS       : "${NODEREQS}
echo "MEMREQS        : "${MEMREQS}
echo -e "\n"

export_vars="\
SAMPLE_YML=$(readlink -e ${SAMPLE_YML}),\
OUT_FOLDER=$(readlink -e ${OUT_FOLDER}),\
RSCRIPTS_DIR=$(readlink -e ${RSCRIPTS_DIR}),\
GENESET_DIR=$(readlink -e ${GENESET_DIR}),\
SPECIES=${SPECIES},\
RERUN_STAGE=${RERUN_STAGE},\
AB_ASSAY_NAME=${AB_ASSAY_NAME}"

UUID=`randomstr 10`

qsub -v ${export_vars} \
     -e ${LOGDIR}/process_10x_seurat_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/process_10x_seurat_${UUID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N process_10x_seurat_${UUID} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
