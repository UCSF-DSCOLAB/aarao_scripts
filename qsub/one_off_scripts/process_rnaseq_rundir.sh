#!/usr/bin/env bash
#
#PBS -V
#PBS -N process_rundir
#PBS -l nodes=1:ppn=1
#PBS -l walltime=3:00:00:00
#PBS -l vmem=20G
#PBS -t 28
#PBS -o /krummellab/data1/arrao/logs/process_rundir.out
#PBS -e /krummellab/data1/arrao/logs/process_rundir.err

set -e
set -o nounset

source /home/arrao/miniconda3/etc/profile.d/conda.sh
conda activate scrnaseq

#python /krummellab/data1/arrao/scripts/python/process_run_dir.py --plate ${PBS_ARRAYID} \
#    --logfile /krummellab/data1/arrao/logs/process_plate${PBS_ARRAYID}_rundir.log \
#    /krummellab/data1/ipi/sequencing/processed/RSQ/plate${PBS_ARRAYID}_rnaseq_new/

python /krummellab/data1/arrao/scripts/python/process_run_dir.py --plate Nova \
    --logfile /krummellab/data1/arrao/logs/process_plateNova_rundir.log \
    /krummellab/data1/ipi/sequencing/processed/RSQ/plateNova_rnaseq_new/

