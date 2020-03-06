#!/bin/bash
set -e
set -o nounset

source /home/arrao/miniconda3/etc/profile.d/conda.sh
conda activate demuxlet

mkdir /scratch/arrao/demuxlet_${SAMPLE} && cd /scratch/arrao/demuxlet_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/demuxlet_${SAMPLE} ; }" EXIT

demuxlet --sam ${SAMFILE} \
         --vcf ${VCF} \
         --field GT \
         --out ${SAMPLE}

mv /scratch/arrao/demuxlet_${SAMPLE} ${OUTDIR}
