#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
INVCF                 : A path to the input vcf
PATIENT               : The name of the patient
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
SNPEFFVERSION         : The version of snpeff to use (v4.3t)
SNPEFFDBDIRECTORY     : A path to the directory storing the snpeff database (/krummellab/data1/ipi/data/databases/snpEff/data)
SNPEFFDB              : The name of the snpEFF database to use (hg38.v85)
SNPEFFDBCONFIG        : A path to a configuration file for the databse (<SNPEFFDBDIRECTORY>/<SNPEFFDB>.config)
NODEREQS              : Default node requirements (default: nodes=1:ppn=1)[overrides global default]
MEMREQS               : Default mem requirements (default: vmem=30gb)[overrides global default]
EOF

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SNPEFFDBDIRECTORY" ]]
then
    SNPEFFDBDIRECTORY="/krummellab/data1/ipi/data/databases/snpEff/data"
else
    SNPEFFDBDIRECTORY=$(readlink -e ${SNPEFFDBDIRECTORY})
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SNPEFFVERSION" ]]
then
    SNPEFFVERSION="v4.3t"
else
    if ! [[ -f /krummellab/data1/ipi/software/snpeff/${SNPEFFVERSION}/snpEff/snpEff.jar ]]
    then
        echo "Could not find jarfile for version ${SNPEFFVERSION} at /krummellab/data1/ipi/software/snpeff/${SNPEFFVERSION}/snpEff/snpEff.jar "
        exit
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SNPEFFDB" ]]
then
    SNPEFFDB="hg38.v85"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "SNPEFFDBCONFIG" ]]
then
    if ! [[ -f ${SNPEFFDBDIRECTORY}/${SNPEFFDB}.config ]]
    then
        echo "Could not find default config file for ${SNPEFFDB} at ${SNPEFFDBDIRECTORY}/${SNPEFFDB}.config "
        exit
    else
        SNPEFFDBCONFIG=${SNPEFFDBDIRECTORY}/${SNPEFFDB}.config
    fi
else
    if ! [[ -f ${SNPEFFDBCONFIG} ]]
    then
        echo "Could not find provided config file ${SNPEFFDBCONFIG} "
        exit
    fi
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=30gb"
fi


echo "Received the following options:"
echo "INVCF                   : "${INVCF-""}
echo "PATIENT                 : "${PATIENT-""}
echo "SNPEFFVERSION           : "${SNPEFFVERSION-""}
echo "SNPEFFDBDIRECTORY       : "${SNPEFFDBDIRECTORY-""}
echo "SNPEFFDB                : "${SNPEFFDB-""}
echo "SNPEFFDBCONFIG          : "${SNPEFFDBCONFIG-""}

if [ ${INVCF-"ERR"} == "ERR" ] || [ ${PATIENT-"ERR"} == "ERR" ] 
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

echo "LOGDIR                  : "${LOGDIR}
echo "NODEREQS                : "${NODEREQS}
echo "MEMREQS                 : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.75 / 1" | bc`

export_vars="\
INVCF=$(readlink -e ${INVCF}),\
PATIENT=${PATIENT},\
SNPEFFVERSION=${SNPEFFVERSION},\
SNPEFFDBDIRECTORY=$(readlink -e ${SNPEFFDBDIRECTORY}),\
SNPEFFDB=${SNPEFFDB},\
SNPEFFDBCONFIG=$(readlink -e ${SNPEFFDBCONFIG}),\
MEMORY=${MEMORY}"

qsub -v ${export_vars} \
     -e ${LOGDIR}/SNPEff_${PATIENT}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/SNPEff_${PATIENT}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N SNPEff_${PATIENT} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
