#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../template_argparse.sh

read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
BAMFILE           : Path to the fastq directory
BARCODELIST       : Path to 10X Barcodes to prioritize
SAMPLE            : The name of the sample
OUTDIR            : A folder within which we will place the output dir
NUMSAMPLES        : The number of samples expected in the pool
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
ONEKGENOMESVCF    : Folder containing per-chromosome 1k genomes vcfs (default: /krummellab/data1/ipi/data/databases/freemuxlet_vcf/ucsc.hg38.liftover.out.nochr.c1_22.nohbb_chromsorted.vcf.gz)
RANDOMSEED        : 21212
NODEREQS          : Default node requirements (default: nodes=1:ppn=10)[overrides global default]
MEMREQS           : Default mem requirements (default: vmem=100gb)[overrides global default]
EOF

if [ ${BAMFILE-"ERR"} == "ERR" ] || [ ${BARCODELIST-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ] || [ ${NUMSAMPLES-"ERR"} == "ERR" ]
then
    echo -e "\nERROR: Required arguments cannot be empty\n"
    print_help  
fi

# Override default NODEREQS and MEMREQS only if the user hasn't specified
if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "NODEREQS" ]]
then
    NODEREQS="nodes=1:ppn=10"
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "MEMREQS" ]]
then
    MEMREQS="vmem=100gb"
fi


if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "ONEKGENOMESVCF" ]]
then
    ONEKGENOMESVCF=/krummellab/data1/ipi/data/databases/freemuxlet_vcf/ucsc.hg38.liftover.out.nochr.c1_22.nohbb_chromsorted.vcf.gz
else
    if [ ! -d ${ONEKGENOMESVCF} ]
    then
       echo -e "\nERROR: ONEKGENOMESVCFDIR dir does not exist: ${ONEKGENOMESVCF}\n" 
    fi
fi

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "RANDOMSEED" ]]
then
    RANDOMSEED=21212
fi

echo "Received the following options:"
echo "BAMFILE            : "${BAMFILE-""}
echo "OUTDIR             : "${OUTDIR-""}
echo "SAMPLE             : "${SAMPLE-""}
echo "NUMSAMPLES         : "${NUMSAMPLES-""}
echo "BARCODELIST        : "${BARCODELIST-""}
echo "ONEKGENOMESVCF     : "${ONEKGENOMESVCF-""}
echo "RANDOMSEED         : "${RANDOMSEED-""}

echo "LOGDIR             : "${LOGDIR}
echo "NODEREQS           : "${NODEREQS}
echo "MEMREQS            : "${MEMREQS}
echo -e "\n"


export_vars="\
BAMFILE=$(readlink -e ${BAMFILE}),\
OUTDIR=$(readlink -f ${OUTDIR}),\
BARCODELIST=${BARCODELIST},\
ONEKGENOMESVCF=$(readlink -e ${ONEKGENOMESVCF}),\
SAMPLE=${SAMPLE},\
NUMSAMPLES=${NUMSAMPLES},\
RANDOMSEED=${RANDOMSEED}"

echo qsub -v ${export_vars} \
     -e ${LOGDIR}/freemuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").err \
     -o ${LOGDIR}/freemuxlet_${SAMPLE}_$(date "+%Y_%m_%d_%H_%M_%S").out \
     -N freemuxlet_${SAMPLE} \
     -l ${NODEREQS} \
     -l ${MEMREQS} \
     $(dirname ${0})/run.sh