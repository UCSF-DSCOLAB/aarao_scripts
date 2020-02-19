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

## Bowtie
### Untrimmed
mkdir -p /scratch/arrao/bowtie_hiseq_untrimmed && cd /scratch/arrao/bowtie_hiseq_untrimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_untrimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_untrimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
    bowtie_hiseq
mv /scratch/arrao/bowtie_hiseq_untrimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/bowtie_hiseq/

mkdir -p /scratch/arrao/bowtie_novaseq_untrimmed && cd /scratch/arrao/bowtie_novaseq_untrimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_untrimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_untrimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
    bowtie_novaseq
mv /scratch/arrao/bowtie_novaseq_untrimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/bowtie_novaseq_untrimmed/
 
### Trimmed
mkdir -p /scratch/arrao/bowtie_hiseq_trimmed && cd /scratch/arrao/bowtie_hiseq_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
    bowtie_hiseq
mv /scratch/arrao/bowtie_hiseq_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/bowtie_hiseq/

mkdir -p /scratch/arrao/bowtie_novaseq_trimmed && cd /scratch/arrao/bowtie_novaseq_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
    bowtie_novaseq
mv /scratch/arrao/bowtie_novaseq_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/bowtie_novaseq_trimmed/

mkdir -p /scratch/arrao/bowtie_novaseq_g_trimmed && cd /scratch/arrao/bowtie_novaseq_g_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --bowtie2 \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_g_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_g_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_bowtie/hg38 \
    bowtie_novaseq
mv /scratch/arrao/bowtie_novaseq_g_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/bowtie_novaseq_g_trimmed/

## STAR
### Untrimmed
mkdir -p /scratch/arrao/star_hiseq_untrimmed && cd /scratch/arrao/star_hiseq_untrimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_untrimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_untrimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_100/hg38 \
    star_hiseq_untrimmed
mv /scratch/arrao/star_hiseq_untrimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/star_hiseq_untrimmed/


mkdir -p /scratch/arrao/star_novaseq_untrimmed && cd /scratch/arrao/star_novaseq_untrimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_untrimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_untrimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_150/hg38 \
    star_novaseq_untrimmed
mv /scratch/arrao/star_novaseq_untrimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/star_novaseq_untrimmed/

### trimmed
mkdir -p /scratch/arrao/star_hiseq_trimmed && cd /scratch/arrao/star_hiseq_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/hiseq_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_100/hg38 \
    star_hiseq_trimmed
mv /scratch/arrao/star_hiseq_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/star_hiseq_trimmed/


mkdir -p /scratch/arrao/star_novaseq_trimmed && cd /scratch/arrao/star_novaseq_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_150/hg38 \
    star_novaseq_trimmed
mv /scratch/arrao/star_novaseq_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/star_novaseq_trimmed/

mkdir -p /scratch/arrao/star_novaseq_g_trimmed && cd /scratch/arrao/star_novaseq_g_trimmed
${RSEM_DIR}/rsem-calculate-expression --paired-end \
    --star \
    --star-gzipped-read-file \
    -p $PBS_NUM_PPN \
    --temporary-folder /scratch/arrao/aligner_test \
    --output-genome-bam \
    --sort-bam-by-coordinate \
    --star-output-genome-bam \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_g_trimmed_1.fastq.gz \
    /data/shared/krummellab/arrao/projects/aligner_comparison/input_fastqs/novaseq_g_trimmed_2.fastq.gz \
    /data/shared/krummellab/ipi/data/refs/rsem/hg38_star_150/hg38 \
    star_novaseq_g_trimmed
mv /scratch/arrao/star_novaseq_g_trimmed/ /data/shared/krummellab/arrao/projects/aligner_comparison/star_novaseq_g_trimmed/


