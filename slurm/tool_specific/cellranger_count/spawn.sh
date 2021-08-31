#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
LIBRARIES_CSV      : Path to a libraries.csv file.
OUTDIR             : A folder within which we will place the output dir
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
CONTAINER          : The singularity container to use (/krummellab/data1/singularity_images/cellranger/6.0.2/cellranger.sif)
SAMPLE             : The name of the sample (parsed from LIBRARIES_CSV by default)
TRANSCRIPTOME      : Path to 10X Indexes (default: /krummellab/data1/ipi/data/refs/10x/refdata-gex-GRCh38-2020-A)
CHEMISTRY          : 10X Chemistry. (default: auto)
FEATUREREF         : Feature reference file for CITE-Seq or hashtags (default: /krummellab/data1/ipi/data/refs/10x/citeseq_resources/TotalSeqC_Human_Universal_Cocktail_v1.csv)
EOF

read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
TIME               : Max Runtime for job (default=2-00:00:00)
CPUSPERTASK        : Default node requirements (default: 12)
MEMPERCPU          : Default mem requirements (default: 13gb)
EOF


if [ ${LIBRARIES_CSV-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if readlink -f ${OUTDIR}
then
    if [[ -d ${OUTDIR} ]] && ! [[ -w ${OUTDIR} ]]
    then
        echo -e "\nERROR: OUTDIR exists but is not writable."
        exit 1
    elif [[ -d $(dirname ${OUTDIR}) ]] && ! [[  -w $(dirname ${OUTDIR}) ]]
    then
        echo -e "\nERROR: Parent dir to OUTDIR exists but is not writable."
        exit 1
    fi
else
    echo -e "\nERROR: parent directory for OUTDIR does not exist\n"
    exit 1
fi



if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SAMPLE" ]]
then
    if grep -q "Gene Expression" ${LIBRARIES_CSV} || [ "$CELLRANGERVERSION" == "3.0.2" ]
    then
        SAMPLE=`grep "Gene Expression" ${LIBRARIES_CSV} | awk -F"," '{print $2}'`
    else
        # Assume the second line has the correct name 
        SAMPLE=`head -2 ${LIBRARIES_CSV} | tail -1 | awk -F"," '{print $2}'`
    fi
fi

# Override defaults only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TIME" ]]
then
    TIME="2-00:00:00"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CPUSPERTASK" ]]
then
    CPUSPERTASK="12"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMPERCPU" ]]
then
    MEMPERCPU="13gb"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TRANSCRIPTOME" ]]
then
    TRANSCRIPTOME=/krummellab/data1/ipi/data/refs/10x/refdata-gex-GRCh38-2020-A
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "FEATUREREF" ]]
then
    FEATUREREF=/krummellab/data1/ipi/data/refs/10x/citeseq_resources/TotalSeqC_Human_Universal_Cocktail_v1.csv
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

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CONTAINER" ]]
then
    CONTAINER="/krummellab/data1/singularity_images/cellranger/6.0.2/cellranger.sif"
else
    if [ ! -f ${CONTAINER} ]
    then
       echo -e "\nERROR: CONTAINER does not exist: ${CONTAINER}\n"
    fi
fi

LOCAL_EXPORT_VARS="\
LIBRARIES_CSV=$(readlink -e ${LIBRARIES_CSV}),\
OUTDIR=$(readlink -f ${OUTDIR}),\
SAMPLE=${SAMPLE},\
CONTAINER=$(readlink -e ${CONTAINER}),\
FEATUREREF=${FEATUREREF},\
TRANSCRIPTOME=$(readlink -e ${TRANSCRIPTOME}),\
CHEMISTRY=${CHEMISTRY}"

JOBNAME=cellranger_count_${SAMPLE}

source $(dirname ${0})/../final_spawn.sh
