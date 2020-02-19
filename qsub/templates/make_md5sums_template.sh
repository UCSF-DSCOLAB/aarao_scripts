#PBS -e <__LOGFOLDER__>/<__PLATE__>_md5sum.err
#PBS -o <__LOGFOLDER__>/<__PLATE__>_md5sum.out
#PBS -N <__PLATE__>_md5sum
#PBS -l <__NODEREQS__>
#PBS -l <__MEMREQS__>

module load CBC python/2.7.3
date
time python /data/shared/krummellab/arrao/scripts/python/process_fastqs.py <__PLATE__>
echo 'Return Status:' $?
