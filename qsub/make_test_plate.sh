#PBS -e /data/shared/krummellab/arrao/logs/make_test_plate.err
#PBS -o /data/shared/krummellab/arrao/logs/make_test_plate.out
#PBS -N make_test_plate
#PBS -l nodes=1:ppn=1
#PBS -l vmem=40gb

#files=($(tree -fi /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq/ | grep R1_001.fastq.gz | sort -r | sed s/"R1_001.fastq.gz"//g))
#for i in {1..50}
#do
#    gunzip -c ${files[i]}R1_001.fastq.gz | head -4000000 | gzip -c > /data/shared/krummellab/ipi/fastqs/PlateTest/test_${i}_R1_001.fastq.gz
#    gunzip -c ${files[i]}R2_001.fastq.gz | head -4000000 | gzip -c > /data/shared/krummellab/ipi/fastqs/PlateTest/test_${i}_R2_001.fastq.gz
#done

files=($(tree -fi /cbc2/data2/samadb/IPI/fastqs/Plate26Rnaseq/ | grep R1_001.fastq.gz | sort -r | sed s/"R1_001.fastq.gz"//g))
for i in {1..50}
do
    j=$(($i + 50))
    gunzip -c ${files[i]}R1_001.fastq.gz | head -4000000 | gzip -c > /data/shared/krummellab/ipi/fastqs/PlateTest/test_${j}_R1_001.fastq.gz
    gunzip -c ${files[i]}R2_001.fastq.gz | head -4000000 | gzip -c > /data/shared/krummellab/ipi/fastqs/PlateTest/test_${j}_R2_001.fastq.gz
done

