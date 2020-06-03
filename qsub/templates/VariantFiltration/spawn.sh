#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
VCFFILE             : A path to the aligned sam/bam/cram file
SAMPLE              : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
INTERVALFILE        : An intervals file to calculate over, use NONE for no itervals (/krummellab/data1/ipi/data/refs/hg38_files/GRCh38_CCDS_sorted.bed.gz)
GATKVERSION         : The version of gatk to use (4.0.2.1)
QUALBYDEPTH         : Filter calls with 
PHREDSCOREDPVAL     : Filter calls with phred scaled p-value [-10log10(p)] greater than this value (23.0)
EOF

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GATKVERSION" ]]
then
    GATKVERSION="4.0.2.1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "INTERVALFILE" ]]
then
    INTERVALFILE="/krummellab/data1/ipi/data/refs/hg38_files/GRCh38_CCDS_sorted.bed.gz"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "PHREDSCOREDPVAL" ]]
then
    PHREDSCOREDPVAL=23.0
fi


echo "Received the following options:"
echo "SAMPLE                  : "${SAMPLE-""}
echo "VCFFILE                 : "${VCFFILE-""}
echo "INTERVALFILE            : "${INTERVALFILE-""}
echo "GATKVERSION             : "${GATKVERSION-""}
echo "PHREDSCOREDPVAL         : "${PHREDSCOREDPVAL-""}

if [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${VCFFILE-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR                  : "${LOGDIR}
echo "NODEREQS                : "${NODEREQS}
echo "MEMREQS                 : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
VCFFILE=$(readlink -f ${VCFFILE}),\
SAMPLE=${SAMPLE},\
INTERVALFILE=${INTERVALFILE},\
MEMORY=${MEMORY},\
GATKVERSION=${GATKVERSION},\
PHREDSCOREDPVAL=${PHREDSCOREDPVAL}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/VariantFiltration_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/VariantFiltration_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N VariantFiltration_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
