#PBS -e /data/shared/krummellab/arrao/logs/star_fusion_darrell.err
#PBS -o /data/shared/krummellab/arrao/logs/star_fusion_darrell.out
#PBS -N star_fusion
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

trap "{ rm -rf /scratch/arrao/* ; }" EXIT

date
source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS

mkdir /scratch/arrao/star-fusion && cd /scratch/arrao/star-fusion

/data/shared/krummellab/ipi/software/STAR-Fusion-v1.5.0/STAR-Fusion --left_fq /cbc2/data2/samadb/IPI/fastqs/Plate12Rnaseq/GC610P12_1/HEP016T1Li_S15_R1_001.fastq.gz \
            --right_fq /cbc2/data2/samadb/IPI/fastqs/Plate12Rnaseq/GC610P12_1/HEP016T1Li_S15_R2_001.fastq.gz \
            --CPU $PBS_NUM_PPN \
            --genome_lib_dir /data/shared/krummellab/ipi/data/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir && \
    mv STAR-Fusion_outdir/ /data/shared/krummellab/arrao/projects/darrell_fusion/HEP016T1Li_fusion/


/data/shared/krummellab/ipi/software/STAR-Fusion-v1.5.0/STAR-Fusion --left_fq /cbc2/data2/samadb/IPI/fastqs/Plate12Rnaseq/GC610P12_3/HEP016T1Tu_S29_R1_001.fastq.gz \
            --right_fq /cbc2/data2/samadb/IPI/fastqs/Plate12Rnaseq/GC610P12_3/HEP016T1Tu_S29_R2_001.fastq.gz \
            --CPU $PBS_NUM_PPN \
            --genome_lib_dir /data/shared/krummellab/ipi/data/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir && \
    mv STAR-Fusion_outdir/ /data/shared/krummellab/arrao/projects/darrell_fusion/HEP016T1Tu_fusion/

echo 'Return Status:' $?