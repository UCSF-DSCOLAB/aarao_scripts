#PBS -e /data/shared/krummellab/arrao/logs/cutadapt.err
#PBS -o /data/shared/krummellab/arrao/logs/cutadapt.out
#PBS -N cutadapt
#PBS -l nodes=1:ppn=24
#PBS -l vmem=40gb

date
/home/shared/cbc/software_cbc/Python-2.7.15/bin/python /data/shared/krummellab/arrao/scripts/python/run_cutadapt_from_yaml.py /data/shared/krummellab/ipi/input_ymls/plates/plate1.rna_paired_align.yml
#/home/arrao/my_usr/bin/tree -fi /data/shared/krummellab/ipi/input_ymls/plates/ | while read yml_file
#do
#  echo "Processing ${yml_file}"
#  time /home/shared/cbc/software_cbc/Python-2.7.15/bin/python run_cutadapt.py ${yml_file}
#done
echo 'Return Status:' $?
