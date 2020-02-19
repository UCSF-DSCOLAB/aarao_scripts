#!/bin/bash
set -e
set -o nounset

source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS

mkdir /scratch/arrao/fastp_${SAMPLE} && cd /scratch/arrao/fastp_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/fastp_${SAMPLE} ; }" EXIT

fastp -i ${FQDIR}/${FQPREFIX}1${FQSUFFIX} \
      -o ${FQPREFIX}1_trimmed${FQSUFFIX} \
      -I ${FQDIR}/${FQPREFIX}2${FQSUFFIX} \
      -O ${FQPREFIX}2_trimmed${FQSUFFIX} \
      --length_required 20 \
      --adapter_sequence CTGTCTCTTATACACATCT \
      --adapter_sequence_r2 CTGTCTCTTATACACATCT \
      --correction  \
      --trim_poly_g  \
      --thread ${PBS_NUM_PPN} \
      -j fastp.json \
      -h fastp.html

mv /scratch/arrao/fastp_${SAMPLE} ${OUTDIR}
