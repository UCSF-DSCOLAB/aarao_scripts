#!/bin/bash
set -e
set -o nounset

mkdir /scratch/arrao/cutadapt_${SAMPLE} && cd /scratch/arrao/cutadapt_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/cutadapt_${SAMPLE}/ ; }" EXIT
            

/home/arrao/.local/bin/cutadapt -j $PBS_NUM_PPN \
                                --trim-n \
                                -a ${READ1ADAPTER} \
                                -A ${READ2ADAPTER} \
                                -o ${SAMPLE}_trimmed_1.fastq.gz \
                                -p ${SAMPLE}_trimmed_2.fastq.gz \
                                -m ${MINREADLEN} \
                                --trim-n ${NEXTSEQTRIM} \
                                ${FQDIR}/${FQPREFIX}1${FQSUFFIX} \
                                ${FQDIR}/${FQPREFIX}2${FQSUFFIX} > log.txt 

mv /scratch/arrao/cutadapt_${SAMPLE} ${OUTDIR}/${SAMPLE}/