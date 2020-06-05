#PBS -e /data/shared/krummellab/arrao/logs/make_ipi_cellranger_indexes.err
#PBS -o /data/shared/krummellab/arrao/logs/make_ipi_cellranger_indexes.out
#PBS -N make_ipi_cellranger_indexes
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

module load CBC cellranger/3.0.0

mkdir /scratch/arrao/cellranger_mkref && cd /scratch/arrao/cellranger_mkref
trap "{ rm -rf /scratch/arrao/cellranger_mkref ; }" EXIT

cellranger mkref --genome=hg38 \
                 --fasta=/cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
                 --genes=/data/shared/krummellab/ipi/data/refs/10x/gtfs/Homo_sapiens.GRCh38.85_10x.gtf \
                 --nthreads=${PBS_NUM_PPN} \
                 --memgb 384

mv hg38/ /data/shared/krummellab/ipi/data/refs/10x/hg38/