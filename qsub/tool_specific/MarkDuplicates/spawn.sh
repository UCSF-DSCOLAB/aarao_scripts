#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
SAMFILE               : A path to the aligned sam/bam/cram file
SAMPLE                : The name of the sample
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
PICARDVERSION         : The version of picard to use (2.22.0)
GENOMEREF             : A path to the reference genome fasta (/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa)
REMOVE_DUPLICATES     : Remove duplicates from the final bam? (false)
VALIDATION_STRINGENCY : Picard validation stringency (default=SILENT)
ASSUME_SORTED         : Assume SAMFILE is sorted? (default=true)
CREATE_INDEX          : Create a bam index? (default=true)
EOF

if [ ${PICARDVERSION-"EMPTY"} == "EMPTY" ]
then
    PICARDVERSION="2.22.0"
fi

if [ ${GENOMEREF-"EMPTY"} == "EMPTY" ]
then
    GENOMEREF="/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa"
fi

if [ ${REMOVE_DUPLICATES-"EMPTY"} == "EMPTY" ]
then
    REMOVE_DUPLICATES="false"
else
    allowed=("true" "false")
    if [[ ! "${allowed[@]}" =~ "${REMOVE_DUPLICATES}" ]]
    then
        echo -e "\nERROR: REMOVE_DUPLICATES must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
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
echo "GENOMEREF               : "${GENOMEREF-""}
echo "REMOVE_DUPLICATES       : "${REMOVE_DUPLICATES-""}
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
SAMFILE=$(readlink -f ${SAMFILE}),
SAMPLE=${SAMPLE},\
GENOMEREF=$(readlink -f ${GENOMEREF}),\
MEMORY=${MEMORY},\
PICARDVERSION=${PICARDVERSION},\
REMOVE_DUPLICATES=${REMOVE_DUPLICATES},
VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY},\
ASSUME_SORTED=${ASSUME_SORTED},\
CREATE_INDEX=${CREATE_INDEX}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/MarkDuplicates_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/MarkDuplicates_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N MarkDuplicates_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
