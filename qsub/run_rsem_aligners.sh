#PBS -e /data/shared/krummellab/arrao/logs/aligner_test.err
#PBS -o /data/shared/krummellab/arrao/logs/aligner_test.out
#PBS -N aligner_test
#PBS -l nodes=1:ppn=48
#PBS -l vmem=384gb

set -o nounset
set -e

date
source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS
RSEM_DIR=/cbc2/data1/lib/rsem/rsem-latest

# mkdir -p /scratch/arrao/test_bowtie_hiseq && cd /scratch/arrao/test_bowtie_hiseq
# ${RSEM_DIR}/rsem-calculate-expression --paired-end \
#     --bowtie2 \
#     -p $PBS_NUM_PPN \
#     --temporary-folder /scratch/arrao/aligner_test \
#     --output-genome-bam \
#     --sort-bam-by-coordinate \
#     /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq/GC0610P22_1/IPIKID100_T1_rna_tumor_S51_R1_001.fastq.gz \
#     /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq/GC0610P22_1/IPIKID100_T1_rna_tumor_S51_R2_001.fastq.gz \
#     /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
#     test_bowtie_hiseq
# mv /scratch/arrao/test_bowtie_hiseq/ /data/shared/krummellab/arrao/projects/aligner_comparison/test_bowtie_hiseq/
# 
# mkdir -p /scratch/arrao/test_bowtie_novaseq && cd /scratch/arrao/test_bowtie_novaseq
# ${RSEM_DIR}/rsem-calculate-expression --paired-end \
#     --bowtie2 \
#     -p $PBS_NUM_PPN \
#     --temporary-folder /scratch/arrao/aligner_test \
#     --output-genome-bam \
#     --sort-bam-by-coordinate \
#     /cbc2/data2/samadb/IPI/fastqs/NovaSeqTrial/GC0610P22_1/IPIKID100_T1_rna_tumor_S99_L004_R1_001.fastq.gz \
#     /cbc2/data2/samadb/IPI/fastqs/NovaSeqTrial/GC0610P22_1/IPIKID100_T1_rna_tumor_S99_L004_R2_001.fastq.gz \
#     /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
#     test_bowtie_novaseq
# mv /scratch/arrao/test_bowtie_novaseq/ /data/shared/krummellab/arrao/projects/aligner_comparison/test_bowtie_novaseq/
 
mkdir -p /scratch/arrao/test_star_hiseq && cd /scratch/arrao/test_star_hiseq
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq/GC0610P22_1/IPIKID100_T1_rna_tumor_S51_R1_001.fastq.gz \
    /cbc2/data2/samadb/IPI/fastqs/Plate22Rnaseq/GC0610P22_1/IPIKID100_T1_rna_tumor_S51_R2_001.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_100/hg38 \
    test_star_hiseq
mv /scratch/arrao/test_star_hiseq/ /data/shared/krummellab/arrao/projects/aligner_comparison/test_star_hiseq/


mkdir -p /scratch/arrao/test_star_novaseq && cd /scratch/arrao/test_star_novaseq
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /cbc2/data2/samadb/IPI/fastqs/NovaSeqTrial/GC0610P22_1/IPIKID100_T1_rna_tumor_S99_L004_R1_001.fastq.gz \
    /cbc2/data2/samadb/IPI/fastqs/NovaSeqTrial/GC0610P22_1/IPIKID100_T1_rna_tumor_S99_L004_R2_001.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_150/hg38 \
    test_star_novaseq
mv /scratch/arrao/test_star_novaseq/ /data/shared/krummellab/arrao/projects/aligner_comparison/test_star_novaseq/


