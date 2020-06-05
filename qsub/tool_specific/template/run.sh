#!/bin/bash
set -e
set -o nounset

module load CBC XXX

mkdir /scratch/arrao/xxx_${SAMPLE} && cd /scratch/arrao/xxx_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/xxx_${SAMPLE} ; }" EXIT


mv /scratch/arrao/xxx_${SAMPLE} ${OUTDIR}/${SAMPLE}/