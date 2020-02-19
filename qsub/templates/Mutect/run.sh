#!/bin/bash
set -e
set -o nounset

module load CBC jdk/7

mkdir /scratch/arrao/MuTect_${SAMPLE} && cd /scratch/arrao/MuTect_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/MuTect_${SAMPLE} ; }" EXIT

out_dir=$(dirname ${BAMFILE})

java -Xmx${MEMORY}g  \
    -Djava.io.tmpdir=/scratch/arrao/${SAMPLE}_javatmp \
     -jar /data/shared/krummellab/ipi/software/muTect-1.1.5/muTect-1.1.5.jar \
     -T MuTect \
     --reference_sequence /data/shared/krummellab/ipi/data/refs/hg38_files/hg38.fa \
     --cosmic /data/shared/krummellab/ipi/data/databases/cosmic/hg38/v88/CosmicCodingMuts.vcf \
     --dbsnp /data/shared/krummellab/ipi/data/databases/dbsnp/hg38/b151/00-All.vcf \
     --input_file:tumor ${BAMFILE} \
     --out mutations.out \
     --vcf mutations.vcf \
     -et NO_ET \
     -K /data/shared/krummellab/ipi/software/muTect-1.1.5/gsamembers_broadinstitute.org.key

mv /scratch/arrao/MuTect_${SAMPLE}/mutations* ${out_dir}/