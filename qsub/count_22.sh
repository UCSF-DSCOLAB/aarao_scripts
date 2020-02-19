#PBS -e /data/shared/krummellab/arrao/logs/count_22.err
#PBS -o /data/shared/krummellab/arrao/logs/count_22.out
#PBS -N count_22
#PBS -l nodes=1:ppn=1
#PBS -l vmem=40gb

date
tree -fi /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq | while read line
do
  echo -e $line `gunzip -c $line|wc -l`
done
echo 'Return Status:' $?
