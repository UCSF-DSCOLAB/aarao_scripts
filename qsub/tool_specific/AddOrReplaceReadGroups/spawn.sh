#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh
source $(dirname ${0})/../../../bash/essentials.sh 

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE               : A path to the aligned sam/bam/cram file
SAMPLE                : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
PICARDVERSION         : The version of picard to use (2.22.0)
RGID                  : Read Group ID (default=<SAMPLE>_1)
RGPL                  : Read Group Platform (default=ILLUMINA)
RGSM                  : Read Group Sample (default=<SAMPLE>)
RGLB                  : Read Group Library (default=<SAMPLE>_1)
RGPU                  : Read Group Platform Unit (default=<random_string>)
VALIDATION_STRINGENCY : Picard validation stringency (default=SILENT)
ASSUME_SORTED         : Assume SAMFILE is sorted? (default=true)
CREATE_INDEX          : Create a bam index? (default=true)
EOF

if [ ${PICARDVERSION-"EMPTY"} == "EMPTY" ]
then
    PICARDVERSION="2.22.0"
fi

if [ ${RGID-"EMPTY"} == "EMPTY" ]
then
    RGID=${SAMPLE}_1
fi

if [ ${RGPL-"EMPTY"} == "EMPTY" ]
then
    RGPL="ILLUMINA"
fi

if [ ${RGSM-"EMPTY"} == "EMPTY" ]
then
    RGSM=${SAMPLE}
fi

if [ ${RGLB-"EMPTY"} == "EMPTY" ]
then
    RGLB=${SAMPLE}_1
fi

if [ ${RGPU-"EMPTY"} == "EMPTY" ]
then
    RGPU=$(randomstr 10)
fi

if [ ${VALIDATION_STRINGENCY-"EMPTY"} == "EMPTY" ]
then
    VALIDATION_STRINGENCY="SILENT"
fi

if [ ${ASSUME_SORTED-"EMPTY"} == "EMPTY" ]
then
    ASSUME_SORTED="true"
else
    allowed=("true" "false")
    if [[ ! "${allowed[@]}" =~ "${ASSUME_SORTED}" ]]
    then
        echo -e "\nERROR: ASSUME_SORTED must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
fi

if [ ${CREATE_INDEX-"EMPTY"} == "EMPTY" ]
then
    CREATE_INDEX="true"
else
    allowed=("true" "false")
    if [[ ! "${allowed[@]}" =~ "${CREATE_INDEX}" ]]
    then
        echo -e "\nERROR: CREATE_INDEX must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
fi


echo "Received the following options:"
echo "SAMPLE                  : "${SAMPLE-""}
echo "SAMFILE                 : "${SAMFILE-""}
echo "PICARDVERSION           : "${PICARDVERSION-""}
echo "RGID                    : "${RGID-""}
echo "RGPL                    : "${RGPL-""}
echo "RGSM                    : "${RGSM-""}
echo "RGLB                    : "${RGLB-""}
echo "RGPU                    : "${RGPU-""}
echo "VALIDATION_STRINGENCY   : "${VALIDATION_STRINGENCY-""}
echo "ASSUME_SORTED           : "${ASSUME_SORTED-""}
echo "CREATE_INDEX            : "${CREATE_INDEX-""}


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
SAMPLE=${SAMPLE},\
SAMFILE=$(readlink -f ${SAMFILE}),SAMPLE=${SAMPLE},\
MEMORY=${MEMORY},\
PICARDVERSION=${PICARDVERSION},\
RGID=${RGID},\
RGPL=${RGPL},\
RGSM=${RGSM},\
RGLB=${RGLB},\
RGPU=${RGPU},\
VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY},\
ASSUME_SORTED=${ASSUME_SORTED},\
CREATE_INDEX=${CREATE_INDEX}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/AddRepRGs_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/AddRepRGs_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N AddRepRGs_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh