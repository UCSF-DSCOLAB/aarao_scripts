#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}


mkdir /scratch/${USER}/SplitNCigarReads_${SAMPLE} && cd /scratch/${USER}/SplitNCigarReads_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/SplitNCigarReads_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})
out_base=$(basename ${SAMFILE%.${extension}})_NCigSplit
gatk --java-options "-Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp -Xmx${MEMORY}g" \
   SplitNCigarReads \
   -R ${GENOMEREF} \
   -I ${SAMFILE} \
   -O ${out_base}.bam

mv /scratch/${USER}/SplitNCigarReads_${SAMPLE}/${out_base}* ${out_dir}/
