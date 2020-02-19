#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

mkdir /scratch/arrao/cellranger_count_${SAMPLEID} && cd /scratch/arrao/cellranger_count_${SAMPLEID} 
trap "{ rm -rf /scratch/arrao/cellranger_count_${SAMPLEID} ; }" EXIT

run_mem=`echo "${MEMORY} * 0.95 / 1" | bc`

cellranger count --id=${SAMPLEID} \
                 --feature-ref /data/shared/krummellab/ipi/data/refs/10x/biolegend_totalseq_hashtags.csv \
                 --libraries /cbc/samadb/ipi_fastqs/10X_3prime/hashed/test1/${SAMPLE}_abc.csv \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${run_mem}

mv /scratch/arrao/cellranger_count_${SAMPLEID}/${SAMPLEID}/ ${OUTDIR}_ABC/

cellranger count --id=${SAMPLEID} \
                 --libraries /cbc/samadb/ipi_fastqs/10X_3prime/hashed/test1/${SAMPLE}_gex.csv \
                 --transcriptome=${TRANSCRIPTOME} \
                 --localcores=${PBS_NUM_PPN} \
                 --localmem=${run_mem}

mv /scratch/arrao/cellranger_count_${SAMPLEID}/${SAMPLEID}/ ${OUTDIR}_GEX/
