#!/bin/bash
set -e
set -o nounset

module load CBC gcc/5.1.0
module load CBC python/2.7.15
module load fastqc/0.11.2

mkdir -p /scratch/arrao/qc_${SAMPLE}/temp && cd /scratch/arrao/qc_${SAMPLE}
trap "{ rm -rf /scratch/arrao/qc_${SAMPLE} ; }" EXIT

# Run tinfoil.cpp first
/data/shared/krummellab/arrao/projects/QC/tinfoil -fq ${FASTQ_FILE} \
                                                 -sample ${SAMPLE} \
                                                 --readlen ${READLEN} \
                                                 --max_records_in_memory 50000000

# Run tinfoil.py
/data/shared/krummellab/arrao/projects/QC/venv/bin/python /data/shared/krummellab/arrao/projects/QC/tinfoil.py ${SAMPLE}_tinfoil/

mv ${SAMPLE}_tinfoil/ ${TINFOILDIR}/${SAMPLE}/

# Now run fastqc
mkdir ${SAMPLE}_fastqc/
fastqc -o ${SAMPLE}_fastqc/ -d temp/ ${FASTQ_FILE}
mv ${SAMPLE}_fastqc/ ${FASTQCDIR}/${SAMPLE}/
