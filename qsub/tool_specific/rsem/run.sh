#!/bin/bash
set -e
set -o nounset

RSEM_DIR=/cbc2/data1/lib/rsem/rsem-latest

mkdir /scratch/arrao/rsem_${SAMPLE} && cd /scratch/arrao/rsem_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/rsem_${SAMPLE} ; }" EXIT

${RSEM_DIR}/rsem-calculate-expression --paired-end \
          -p ${PBS_NUM_PPN} \
          --alignments ${BAMFILE} \
          --no-bam-output \
          --temporary-folder /scratch/arrao/rsem_${SAMPLE}/temp \
          ${RSEMREFERENCE} \
          ${SAMPLE}

mv /scratch/arrao/rsem_${SAMPLE} ${OUTDIR}/