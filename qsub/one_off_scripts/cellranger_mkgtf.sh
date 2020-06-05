#PBS -e /data/shared/krummellab/arrao/logs/make_ipi_cellranger_indexes.err
#PBS -o /data/shared/krummellab/arrao/logs/make_ipi_cellranger_indexes.out
#PBS -N make_ipi_cellranger_indexes
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

module load CBC cellranger/3.0.0

cellranger mkgtf /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf /data/shared/krummellab/ipi/data/refs/10x/hg38/Homo_sapiens.GRCh38.85_10x.gtf