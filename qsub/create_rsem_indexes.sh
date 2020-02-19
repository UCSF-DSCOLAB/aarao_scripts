#PBS -e /data/shared/krummellab/arrao/logs/star_index.err
#PBS -o /data/shared/krummellab/arrao/logs/star_index.out
#PBS -N star_index
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

date

source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS
RSEM_DIR=/cbc2/data1/lib/rsem/rsem-latest

${RSEM_DIR}/rsem-prepare-reference --gtf /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38


${RSEM_DIR}/rsem-prepare-reference --gtf /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
    --star \
    --star-sjdboverhang 100 \
    -p $PBS_NUM_PPN \
    /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_100/hg38


${RSEM_DIR}/rsem-prepare-reference --gtf /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
    --star \
    --star-sjdboverhang 150 \
    -p $PBS_NUM_PPN \
    /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_150/hg38