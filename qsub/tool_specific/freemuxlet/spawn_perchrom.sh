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
EOF

read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
ONEKGENOMESVCFDIR : Folder containing per-chromosome 1k genomes vcfs (default: /krummellab/data1/ipi/data/databases/1k_genomes/genic_chromosomal_vcfs)
NODEREQS          : Default node requirements (default: nodes=1:ppn=10)[overrides global default]
MEMREQS           : Default mem requirements (default: vmem=100gb)[overrides global default]
EOF

if [ ${BAMFILE-"ERR"} == "ERR" ] || [ ${BARCODELIST-"ERR"} == "ERR" ] || [ ${SAMPLE-"ERR"} == "ERR" ] || [ ${OUTDIR-"ERR"} == "ERR" ]
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


if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "ONEKGENOMESVCFDIR" ]]
then
    ONEKGENOMESVCFDIR=/krummellab/data1/ipi/data/databases/1k_genomes/genic_chromosomal_vcfs
else
    if [ ! -d ${ONEKGENOMESVCFDIR} ]
    then
       echo -e "\nERROR: ONEKGENOMESVCFDIR dir does not exist: ${ONEKGENOMESVCFDIR}\n" 
    fi
fi

if [ -d ${OUTDIR} ]
then
    echo "OUTDIR exists. Please choose an OUTDIR that does not exist."
    exit
else
    mkdir -p ${OUTDIR}
fi

echo "Received the following options:"
echo "BAMFILE            : "${BAMFILE-""}
echo "OUTDIR             : "${OUTDIR-""}
echo "SAMPLE             : "${SAMPLE-""}
echo "BARCODELIST        : "${BARCODELIST-""}
echo "ONEKGENOMESVCFDIR  : "${ONEKGENOMESVCFDIR-""}

echo "LOGDIR             : "${LOGDIR}
echo "NODEREQS           : "${NODEREQS}
echo "MEMREQS            : "${MEMREQS}
echo -e "\n"

MEMORY=`echo "$(echo ${MEMREQS} | sed 's/[^0-9]*//g')*0.9 / 1" | bc`


export_vars="\
BAMFILE=$(readlink -e ${BAMFILE}),\
OUTDIR=$(readlink -e ${OUTDIR}),\
BARCODELIST=${BARCODELIST}"

for chrom in {1..22} X Y MT
do
    chrom_export_vars=",\
ONEKGENOMESVCF=$(readlink -e ${ONEKGENOMESVCFDIR})/${chrom}_genic.vcf.gz,\
SAMPLE=${SAMPLE}_${chrom}"

    qsub -v ${export_vars}${chrom_export_vars} \
         -e ${LOGDIR}/freemuxlet_bychrom_${SAMPLE}_${chrom}_$(date "+%Y_%m_%d_%H_%M_%S").err \
         -o ${LOGDIR}/freemuxlet_bychrom_${SAMPLE}_${chrom}_$(date "+%Y_%m_%d_%H_%M_%S").out \
         -N freemuxlet_bychrom_${SAMPLE}_${chrom} \
         -l ${NODEREQS} \
         -l ${MEMREQS} \
         $(dirname ${0})/run.sh
    break
done