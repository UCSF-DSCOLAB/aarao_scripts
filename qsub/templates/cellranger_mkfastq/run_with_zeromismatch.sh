#!/bin/bash
set -e
set -o nounset

module load CBC cellranger/3.0.2

scratch_dir=/scratch/arrao/cellranger_mkfastq_${UUID}
mkdir -p ${scratch_dir}/working && cd ${scratch_dir}/working
trap "{ rm -rf ${scratch_dir} ; }" EXIT

run_mem=`echo "${MEMORY} * 0.95 / 1" | bc`

cellranger mkfastq --csv=${SAMPLESHEET} \
                   --run=${BCLDIR} \
                   --barcode-mismatches=0 \
                   --lanes=1,2 \
                   --output-dir=${scratch_dir}/fastqs \
                   --localcores=${PBS_NUM_PPN} \
                   --localmem=${run_mem}

find ${scratch_dir}/fastqs/ -name "*.gz" -exec dirname {} \; | sort -u | while read subfolder
do
  cd ${subfolder}
  echo "Validating entries within ${subfolder}"
  echo "Testing Gzip integrity"
  ls *gz | xargs -n1 /krummellab/data1/arrao/scripts/bash/on_path/test_pigz_integrity
  echo "Generating md5sums"
  md5sum *.fastq.gz > md5checksum.txt
  echo ""
done

echo "Filesizes for All generated fastqs:"
find ${scratch_dir}/fastqs -name "*.gz" -exec du -h {} \; | sort -k2d

mv ${scratch_dir}/fastqs ${OUTDIR}/
