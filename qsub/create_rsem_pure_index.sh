#PBS -e /data/shared/krummellab/arrao/logs/create_pure_rsem_index.err
#PBS -o /data/shared/krummellab/arrao/logs/create_pure_rsem_index.out
#PBS -N create_pure_rsem_index
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

RSEM_DIR=/cbc2/data1/lib/rsem/rsem-latest

mkdir /data/shared/krummellab/ipi/data/refs/rsem/hg38_pure_rsem/
${RSEM_DIR}/rsem-prepare-reference --gtf /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
    -p ${PBS_NUM_PPN} \
    /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_pure_rsem/hg38
