#!/bin/bash
source /krummellab/data1/ipi/software/samtools/samtools-1.10-usr/SOURCE_THIS
source $(dirname ${0})/../template_argparse.sh
set -e
set -o nounset

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
TUMOR_BAMFILE       : A path to the aligned tumor bam file
NORMAL_BAMFILE      : A path to the aligned normal bam file
PATIENT             : The name of the sample
OUTDIR              : A path to a folder where outputs will be written. Will be created if requried
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
TUMOR_SAMPLE_NAME   : The tumor sample name to pull from the bamfile. (the first SN pulled from <TUMOR_BAMFILE>)
NORMAL_SAMPLE_NAME  : The normal sample name to pull from the bamfile. (the first SN pulled from <NORMAL_BAMFILE>)
GATKVERSION         : The version of gatk to use (4.0.2.1)
GENOMEREF           : A path to the reference genome fasta (/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa)
TARGET_INTERVALS    : A path to a file containing intervals to search over (/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list)
GERMLINE_RESOURCE   : The germline resource to use as background (/krummellab/data1/ipi/data/databases/gnomAD/ExAC/ExAC.r1.sites.hg38.liftoverFromhg19.vep.chromsorted.vcf.gz)
NODEREQS            : Default node requirements (default: nodes=1:ppn=4)[overrides global default]
MEMREQS             : Default mem requirements (default: vmem=75gb)[overrides global default]
EOF

echo "Received the following options:"
echo "TUMOR_BAMFILE       : "${TUMOR_BAMFILE-""}
echo "NORMAL_BAMFILE      : "${NORMAL_BAMFILE-""}
echo "PATIENT             : "${PATIENT-""}
echo "OUTDIR              : "${OUTDIR-""}

if [ ${TUMOR_BAMFILE-"ERR"} == "ERR" ] || [ ${NORMAL_BAMFILE-"ERR"} == "ERR" ] || [ ${PATIENT-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GATKVERSION" ]]
then
    GATKVERSION="4.0.2.1"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GENOMEREF" ]]
then
    GENOMEREF="/krummellab/data1/ipi/data/refs/hg38_files/hg38.fa"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TARGET_INTERVALS" ]]
then
    TARGET_INTERVALS="/krummellab/data1/ipi/data/refs/AgilentSureSelect_Human_All_Exome_V6/hg38/S07604514_Regions_nochr.interval_list"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "GERMLINE_RESOURCE" ]]
then
    GERMLINE_RESOURCE="/krummellab/data1/ipi/data/databases/gnomAD/ExAC/ExAC.r1.sites.hg38.liftoverFromhg19.vep.chromsorted.vcf.gz"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "TUMOR_SAMPLE_NAME" ]]
then
    temp=$(samtools view -H --no-PG ${TUMOR_BAMFILE} | grep -m1 "@RG")
    temp=($(echo ${temp#*SM:}))
    TUMOR_SAMPLE_NAME=${temp[0]}
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NORMAL_SAMPLE_NAME" ]]
then
    temp=$(samtools view -H --no-PG ${NORMAL_BAMFILE} | grep -m1 "@RG")
    temp=($(echo ${temp#*SM:}))
    NORMAL_SAMPLE_NAME=${temp[0]}
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=4"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=75gb"
fi


echo "TUMOR_SAMPLE_NAME   : "${TUMOR_SAMPLE_NAME-""}
echo "NORMAL_SAMPLE_NAME  : "${NORMAL_SAMPLE_NAME-""}
echo "GATKVERSION         : "${GATKVERSION-""}
echo "GENOMEREF           : "${GENOMEREF-""}
echo "TARGET_INTERVALS    : "${TARGET_INTERVALS-""}
echo "GERMLINE_RESOURCE   : "${GERMLINE_RESOURCE-""}
echo "LOGDIR              : "${LOGDIR}
echo "NODEREQS            : "${NODEREQS}
echo "MEMREQS             : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`

export_vars="\
TUMOR_BAMFILE=$(readlink -f ${TUMOR_BAMFILE}),\
NORMAL_BAMFILE=$(readlink -f ${NORMAL_BAMFILE}),\
PATIENT=${PATIENT},\
OUTDIR=${OUTDIR},\
TUMOR_SAMPLE_NAME=${TUMOR_SAMPLE_NAME},\
NORMAL_SAMPLE_NAME=${NORMAL_SAMPLE_NAME},\
GATKVERSION=${GATKVERSION},\
GENOMEREF=$(readlink -f ${GENOMEREF}),\
TARGET_INTERVALS=$(readlink -f ${TARGET_INTERVALS}),\
GERMLINE_RESOURCE=$(readlink -f ${GERMLINE_RESOURCE}),\
MEMORY=${MEMORY}"

qsub -v  ${export_vars} \
     -e ${LOGDIR}/MuTect2_${PATIENT}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/MuTect2_${PATIENT}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N MuTect2_${PATIENT} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh
