#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
LIBRARIES_CSV  : Path to a libraries.csv file.
OUTDIR         : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
SAMPLE         : The name of the sample (parsed from LIBRARIES_CSV by default)
NODEREQS       : Default node requirements (default: nodes=1:ppn=32)[overrides global default]
MEMREQS        : Default mem requirements (default: vmem=250gb)[overrides global default]
TRANSCRIPTOME  : Path to 10X Indexes (default: /krummellab/data1/ipi/data/refs/10x/hg38)
CHEMISTRY      : 10X Chemistry. (default: auto)
FEATUREREF     : Feature reference file for CITE-Seq or hashtags (default: /krummellab/data1/ipi/data/refs/10x/biolegend_totalseq_hashtags.csv)
EOF

if [ ${LIBRARIES_CSV-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi


if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SAMPLE" ]]
then
    SAMPLE=`grep "Gene Expression" ${LIBRARIES_CSV} | awk -F"," '{print $2}'`
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

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TRANSCRIPTOME" ]]
then
    TRANSCRIPTOME=/krummellab/data1/ipi/data/refs/10x/hg38
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "FEATUREREF" ]]
then
    FEATUREREF=/krummellab/data1/ipi/data/refs/10x/biolegend_totalseq_hashtags.csv
else
    if [ ! -f ${FEATUREREF} ]
    then
       echo -e "\nERROR: FEATUREREF file does not exist: ${FEATUREREF}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CHEMISTRY" ]]
then
    CHEMISTRY="auto"
else
    allowed=('auto' 'threeprime' 'fiveprime' 'SC3Pv1' 'SC3Pv2' 'SC3Pv3' 'SC5P-PE' 'SC5P-R2')
    if [[ ! "${allowed[@]}" =~ "${CHEMISTRY}" ]]
    then
        echo -e "\nERROR: CHEMISTRY must be one of \n\t${allowed[@]}\n"
        exit 1
    fi
fi 


echo "Received the following options:"
echo "LIBRARIES_CSV  : "${LIBRARIES_CSV-""}
echo "OUTDIR         : "${OUTDIR-""}
echo "SAMPLE         : "${SAMPLE-""}
echo "TRANSCRIPTOME  : "${TRANSCRIPTOME-""}
echo "FEATUREREF     : "${FEATUREREF-""}
echo "CHEMISTRY      : "${CHEMISTRY-""}

echo "LOGDIR   : "${LOGDIR}
echo "NODEREQS : "${NODEREQS}
echo "MEMREQS  : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
LIBRARIES_CSV=$(readlink -e ${LIBRARIES_CSV}),\
SAMPLE=${SAMPLE},\
OUTDIR=${OUTDIR},\
FEATUREREF=${FEATUREREF},\
TRANSCRIPTOME=$(readlink -e ${TRANSCRIPTOME}),\
CHEMISTRY=${CHEMISTRY},\
MEMORY=${MEMORY}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/cellranger_count_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/cellranger_count_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N cellranger_count_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
