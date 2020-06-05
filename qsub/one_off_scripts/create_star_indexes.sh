#PBS -e /data/shared/krummellab/arrao/logs/create_star_index.err
#PBS -o /data/shared/krummellab/arrao/logs/create_star_index.out
#PBS -N star_index
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS

mkdir /scratch/arrao/star_index && cd /scratch/arrao/star_index
trap "{ rm -rf /scratch/arrao/star_index ; }" EXIT

mkdir hg38_100sjdb/
STAR --runThreadN $PBS_NUM_PPN \
     --runMode genomeGenerate \
     --genomeDir hg38_100sjdb/ \
     --genomeFastaFiles /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
     --sjdbGTFfile /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
     --sjdbOverhang 100
mv hg38_100sjdb/ /data/shared/krummellab/ipi/data/refs/star/hg38_100sjdb/

mkdir hg38_150sjdb/
STAR --runThreadN $PBS_NUM_PPN \
     --runMode genomeGenerate \
     --genomeDir hg38_150sjdb/ \
     --genomeFastaFiles /cbc/resources/human/hg38/bwa_index/Homo_sapiens.GRCh38.dna.wholechrom.fa \
     --sjdbGTFfile /cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf \
     --sjdbOverhang 150
mv hg38_150sjdb/ /data/shared/krummellab/ipi/data/refs/star/hg38_150sjdb/