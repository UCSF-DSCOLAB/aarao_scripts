#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/arrao/scripts/bash/essentials.sh
uuid=`randomstr 10`

module load CBC cellranger/3.0.2

mkdir -p /scratch/arrao/cellranger_mkfastq_${SAMPLE}_${uuid}/working  /scratch/arrao/cellranger_mkfastq_${SAMPLE}_${uuid}/working && cd /scratch/arrao/cellranger_mkfastq_${SAMPLE}_${uuid}/working 
trap "{ rm -rf /scratch/arrao/cellranger_mkfastq_${SAMPLE}_${uuid} ; }" EXIT

fastq_dir=$(dirname $(pwd))/fastqs

lanes_argstring=" "
if [ ${LANES} != "ALL" ]
then
    lanes_argstring="--lanes=${LANES} "
fi

cellranger mkfastq --csv=${SAMPLESHEET} \
                   --run=${BCLDIR} \
                   --barcode-mismatches=0 \
                   ${lanes_argstring} \
                   --output-dir=${fastq_dir} \
                   --localcores=${PBS_NUM_PPN} \
                   --localmem=${MEMORY}


find ${fastq_dir}/ -name "*.gz" -exec dirname {} \; | sort -u | while read subfolder
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
find ${fastq_dir}/fastqs -name "*.gz" -exec du -h {} \; | sort -k2d


if [ -f ${OUTDIR} ]
then
    mv ${fastq_dir} ${OUTDIR}/${SAMPLE}_mkfastq_${uuid}
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${fastq_dir} ${OUTDIR}
fi
