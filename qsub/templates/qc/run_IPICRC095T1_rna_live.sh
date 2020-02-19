#!/bin/bash
#PBS -e /data/shared/krummellab/ipi/qc/tinfoil/Plate25Rnaseq/logs/tinfoil_IPICRC095T1_rna_live_S8_R12_001.err
#PBS -o /data/shared/krummellab/ipi/qc/tinfoil/Plate25Rnaseq/logs/tinfoil_IPICRC095T1_rna_live_S8_R12_001.out
#PBS -N tinfoil_25
#PBS -l nodes=1:ppn=1
#PBS -l vmem=200gb

set -e
set -o nounset

module load CBC gcc/5.1.0
module load CBC python/2.7.15


mkdir /scratch/arrao/tinfoil_XXX && cd /scratch/arrao/tinfoil_XXX
trap "{ rm -rf /scratch/arrao/tinfoil_XXX ; }" EXIT

/data/shared/krummellab/arrao/projects/QC/tinfoil -fq /cbc2/data2/samadb/IPI/fastqs/Plate25Rnaseq/GC0610P25_3/IPICRC095T1_rna_live_S8_R1_001.fastq.gz \
                                                 -sample IPICRC095T1_rna_live_S8_R1_001 \
                                                 --readlen 101 \
                                                 --max_records_in_memory 100000000
/data/shared/krummellab/arrao/projects/QC/venv/bin/python /data/shared/krummellab/arrao/projects/QC/tinfoil.py IPICRC095T1_rna_live_S8_R1_001_tinfoil/
mv IPICRC095T1_rna_live_S8_R1_001_tinfoil/ /data/shared/krummellab/ipi/qc/tinfoil/Plate25Rnaseq/IPICRC095T1_rna_live_S8_R1_001/


/data/shared/krummellab/arrao/projects/QC/tinfoil -fq /cbc2/data2/samadb/IPI/fastqs/Plate25Rnaseq/GC0610P25_3/IPICRC095T1_rna_live_S8_R2_001.fastq.gz \
                                                 -sample IPICRC095T1_rna_live_S8_R2_001 \
                                                 --readlen 101 \
                                                 --max_records_in_memory 100000000
/data/shared/krummellab/arrao/projects/QC/venv/bin/python /data/shared/krummellab/arrao/projects/QC/tinfoil.py IPICRC095T1_rna_live_S8_R2_001_tinfoil/
mv IPICRC095T1_rna_live_S8_R2_001_tinfoil/ /data/shared/krummellab/ipi/qc/tinfoil/Plate25Rnaseq/IPICRC095T1_rna_live_S8_R2_001/
