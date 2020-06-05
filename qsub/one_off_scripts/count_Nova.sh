#PBS -e /data/shared/krummellab/arrao/logs/count_Nova.err
#PBS -o /data/shared/krummellab/arrao/logs/count_Nova.out
#PBS -N count_Nova
#PBS -l nodes=1:ppn=1
#PBS -l vmem=40gb

date
tree -fi /cbc2/data2/samadb/IPI/fastqs/PlateNovaRnaseq | while read line
do
  echo -e $line `gunzip -c $line|wc -l`
done
echo 'Return Status:' $?
