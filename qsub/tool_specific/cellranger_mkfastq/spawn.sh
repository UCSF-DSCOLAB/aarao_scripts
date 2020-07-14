#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
BCLDIR             : path to a directory containing BCL files
OUTDIR             : A folder within which we will place the output dir
FLOWCELLID         : A unique ID for the flowcell being processed.
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
CELLRANGERVERSION  : The version of cellranger to use (3.0.2)
BARCODEMISMATCHES  : Number of mismatches to allow in the barcode. (default=1)
LANES              : Comma separated list of lanes to use. NO SPACES. (default=ALL)
SAMPLESHEET        : Path to an Illumina samplesheet (default=<BCLDIR>/SampleSheet.csv)
NODEREQS           : Default node requirements (default: nodes=1:ppn=32)[overrides global default]
MEMREQS            : Default mem requirements (default: vmem=250gb)[overrides global default]
EOF

if [ ${BCLDIR-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ] || [ ${FLOWCELLID-"ERR"} == "ERR" ]
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

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SAMPLESHEET" ]]
then
    SAMPLESHEET=${BCLDIR}/SampleSheet.csv
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "BARCODEMISMATCHES" ]]
then
    BARCODEMISMATCHES=1
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "LANES" ]]
then
    LANES=ALL
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

echo "Received the following options:"
echo "BCLDIR              : "${BCLDIR-""}
echo "OUTDIR              : "${OUTDIR-""}
echo "FLOWCELLID          : "${FLOWCELLID-""}

echo "CELLRANGERVERSION   : "${CELLRANGERVERSION-""}
echo "SAMPLESHEET         : "${SAMPLESHEET-""}
echo "BARCODEMISMATCHES   : "${BARCODEMISMATCHES-""}
echo "LANES               : "${LANES-""}
echo "LOGDIR              : "${LOGDIR}
echo "NODEREQS            : "${NODEREQS}
echo "MEMREQS             : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
BCLDIR=$(readlink -e ${BCLDIR}),\
SAMPLESHEET=$(readlink -e ${SAMPLESHEET}),\
OUTDIR=${OUTDIR},\
FLOWCELLID=${FLOWCELLID},\
CELLRANGERVERSION=${CELLRANGERVERSION},\
BARCODEMISMATCHES=${BARCODEMISMATCHES},\
LANES=${LANES},\
MEMORY=${MEMORY}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/cellranger_mkfastq_${FLOWCELLID}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/cellranger_mkfastq_${FLOWCELLID}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_mkfastq_${FLOWCELLID} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
