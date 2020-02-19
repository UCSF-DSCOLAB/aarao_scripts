#PBS -e /data/shared/krummellab/arrao/logs/bwa_index_rals.err
#PBS -o /data/shared/krummellab/arrao/logs/bwa_index_rals.out
#PBS -N bwa_index_rals
#PBS -l nodes=1:ppn=1
#PBS -l vmem=100gb

module load CBC bwa
date
bwa index /data/shared/krummellab/ipi/data/refs/bwa/human.rrna.fa
echo 'Return Status:' $?
