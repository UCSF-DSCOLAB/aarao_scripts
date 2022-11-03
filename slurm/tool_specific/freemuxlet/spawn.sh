#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
BAMFILE                  : Path to the fastq directory
BARCODELIST              : Path to 10X Barcodes to prioritize
SAMPLE                   : The name of the sample
OUTDIR                   : A folder within which we will place the output dir
NUMSAMPLES               : The number of samples expected in the pool
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
CONTAINER                : The singularity container to use (/krummellab/data1/singularity_images/popscle/da70fc78da385ef049e0e890342acfd62842cae0/popscle.sif)
POPSCLE_HELPER_TOOLS_DIR : A path to the popscle_helper_tools folder (/krummellab/data1/ipi/software/popscle_helper_tools)
ONEKGENOMESVCF           : Folder containing per-chromosome 1k genomes vcfs (default: /krummellab/data1/ipi/data/databases/freemuxlet_vcf/ucsc.hg38.liftover.out.withchr.c1_22.nohbb.vcf.gz)
RANDOMSEED               : 21212
NO_TAG_UMI               : Positional arg to skip the \`--tag-UMI\` option.
DSC_ONLY                 : Positional arg to run only DSC and not FMX, used pre-merge
FMX_ONLY                 : Positional arg to run only FMX and not DSC before, used post merge
EOF

read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
CPUSPERTASK              : Default node requirements (default: 10)
MEMPERCPU                : Default mem requirements (default: 10gb)
EOF

if [ ${BAMFILE-"ERR"} == "ERR" ] || [ ${BARCODELIST-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ] || [ ${NUMSAMPLES-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

# Override defaults only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CPUSPERTASK" ]]
then
    CPUSPERTASK="10"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMPERCPU" ]]
then
    MEMPERCPU="10gb"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "CONTAINER" ]]
then
    CONTAINER="/krummellab/data1/singularity_images/popscle/da70fc78da385ef049e0e890342acfd62842cae0/popscle.sif"
else
    if [ ! -f ${CONTAINER} ]
    then
       echo -e "\nERROR: CONTAINER does not exist: ${CONTAINER}\n"
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "POPSCLE_HELPER_TOOLS_DIR" ]]
then
    POPSCLE_HELPER_TOOLS_DIR=/krummellab/data1/ipi/software/popscle_helper_tools
else
    if [ ! -d ${POPSCLE_HELPER_TOOLS_DIR} ]
    then
       echo -e "\nERROR: POPSCLE_HELPER_TOOLS_DIR dir does not exist: ${POPSCLE_HELPER_TOOLS_DIR}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "ONEKGENOMESVCF" ]]
then
    ONEKGENOMESVCF=/krummellab/data1/ipi/data/databases/freemuxlet_vcf/ucsc.hg38.liftover.out.withchr.c1_22.nohbb.vcf.gz
else
    if [ ! -f ${ONEKGENOMESVCF} ]
    then
       echo -e "\nERROR: ONEKGENOMESVCFDIR dir does not exist: ${ONEKGENOMESVCF}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "RANDOMSEED" ]]
then
    RANDOMSEED=21212
fi

if [[ "${POSITIONAL_ARGS[@]}" =~ "NO_TAG_UMI" ]]
then
    NO_TAG_UMI=TRUE
else
    NO_TAG_UMI=FALSE
fi


if [[ "${POSITIONAL_ARGS[@]}" =~ "DSC_ONLY" ]]
then
    DSC_ONLY=TRUE
else
    DSC_ONLY=FALSE
fi

if [[ "${POSITIONAL_ARGS[@]}" =~ "FMX_ONLY" ]]
then
    FMX_ONLY=TRUE
else
    FMX_ONLY=FALSE
fi


LOCAL_EXPORT_VARS="\
BAMFILE=$(readlink -e ${BAMFILE}),\
OUTDIR=$(readlink -f ${OUTDIR}),\
BARCODELIST=${BARCODELIST},\
CONTAINER=$(readlink -e ${CONTAINER}),\
POPSCLE_HELPER_TOOLS_DIR=$(readlink -e ${POPSCLE_HELPER_TOOLS_DIR}),\
ONEKGENOMESVCF=$(readlink -e ${ONEKGENOMESVCF}),\
SAMPLE=${SAMPLE},\
NUMSAMPLES=${NUMSAMPLES},\
RANDOMSEED=${RANDOMSEED},\
NO_TAG_UMI=${NO_TAG_UMI},\
DSC_ONLY=${DSC_ONLY},\
FMX_ONLY=${FMX_ONLY}"

JOBNAME=freemuxlet_${SAMPLE}

source $(dirname ${0})/../final_spawn.sh
