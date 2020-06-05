#!/bin/bash
set -e
set -o nounset

source /home/arrao/miniconda3/etc/profile.d/conda.sh
conda activate scrnaseq

mkdir /scratch/arrao/citeseq_count_${SAMPLE} && cd /scratch/arrao/citeseq_count_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/citeseq_count_${SAMPLE} ; }" EXIT

if [ ${CHEM} == "V2" ]
then
    UMIL=26
else
    UMIL=28
fi

CELLS=`wc -l ${WHITELIST} | awk '{print $1}'`

CITE-seq-Count -T ${PBS_NUM_PPN} \
               -R1 ${FQ1} \
               -R2 ${FQ2} \
               -t  ${CSTAGS} \
               -cbf 1 -cbl 16 \
               -umif 17 -umil ${UMIL} \
               -cells ${CELLS} \
               -wl ${WHITELIST}
               # -n 400

mv Results/ ${OUTDIR}/
