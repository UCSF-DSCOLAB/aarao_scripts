#PBS -e /data/shared/krummellab/arrao/logs/27_md5sum.err
#PBS -o /data/shared/krummellab/arrao/logs/27_md5sum.out
#PBS -N 27_md5sum
#PBS -l nodes=1:ppn=24
#PBS -l vmem=20gb

module load CBC python/2.7.3
date
time python /data/shared/krummellab/arrao/scripts/python/process_fastqs.py --fastq_dir /data/shared/krummellab/ipi/fastqs 27
echo 'Return Status:' $?

