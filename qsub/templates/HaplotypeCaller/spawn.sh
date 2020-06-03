#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE             : A path to the aligned sam/bam/cram file
SAMPLE              : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
INTERVALFILE        : An intervals file to calculate over, use NONE for no itervals (/krummellab/data1/ipi/data/refs/hg38_files/GRCh38_CCDS_sorted.bed.gz)
GATKVERSION         : The version of gatk to use (4.0.2.1)
GENOMEREF           : A path to the reference genome fasta (/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa)
GVCF                : Use GVCF filtering? (false)
STANDCALLCONF       : The minimum phred-scaled confidence threshold at which variants should be called (20.0)
USESOFTCLIPPEDBASES : Analyze soft clipped bases in the reads? (false)
EOF

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GATKVERSION" ]]
then
    GATKVERSION="4.0.2.1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GENOMEREF" ]]
then
    GENOMEREF="/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "INTERVALFILE" ]]
then
    INTERVALFILE="/krummellab/data1/ipi/data/refs/hg38_files/GRCh38_CCDS_sorted.bed.gz"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "DBSNPFILE" ]]
then
    DBSNPFILE="/krummellab/data1/ipi/data/databases/dbsnp/hg38/b151/00-All.vcf.gz"
fi

if [[ ! "${POSITIONAL_ARGS[@]}" =~ "GVCF" ]]
then
    GVCF="FALSE"
else
    GVCF="TRUE"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "STANDCALLCONF" ]]
then
    STANDCALLCONF=20.0
else
    if (( $(echo "$STANDCALLCONF <= 0.00 || 60.0 <= $STANDCALLCONF" |bc -l) ))
    then
        echo -e "\nERROR: STANDCALLCONF must be in the range (0.00, 60.0) \n"
        print_help
    fi
fi

if [[ ! "${POSITIONAL_ARGS[@]}" =~ "USESOFTCLIPPEDBASES" ]]
then
    DONTUSESOFTCLIPPEDBASES="true"
else
    DONTUSESOFTCLIPPEDBASES="false"
fi


echo "Received the following options:"
echo "SAMPLE                  : "${SAMPLE-""}
echo "SAMFILE                 : "${SAMFILE-""}
echo "INTERVALFILE            : "${INTERVALFILE-""}
echo "DBSNPFILE               : "${DBSNPFILE-""}
echo "GATKVERSION             : "${GATKVERSION-""}
echo "GENOMEREF               : "${GENOMEREF-""}
echo "GVCF                    : "${GVCF-""}
echo "STANDCALLCONF           : "${STANDCALLCONF-""}
echo "DONTUSESOFTCLIPPEDBASES : "${DONTUSESOFTCLIPPEDBASES-""}

if [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${SAMFILE-"ERR"} == "ERR" ]
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
SAMFILE=$(readlink -f ${SAMFILE}),\
SAMPLE=${SAMPLE},\
GENOMEREF=$(readlink -f ${GENOMEREF}),\
INTERVALFILE=${INTERVALFILE},\
DBSNPFILE=${DBSNPFILE},\
MEMORY=${MEMORY},\
GATKVERSION=${GATKVERSION},\
GVCF=${GVCF},\
STANDCALLCONF=${STANDCALLCONF},\
DONTUSESOFTCLIPPEDBASES=${DONTUSESOFTCLIPPEDBASES}"


qsub -v ${export_vars} \
     -e ${LOGDIR}/HaplotypeCaller_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/HaplotypeCaller_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N HaplotypeCaller_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
