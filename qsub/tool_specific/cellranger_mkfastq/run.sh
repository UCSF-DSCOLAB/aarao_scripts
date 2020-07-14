#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/${USER}/aarao_scripts/bash/essentials.sh
uuid=`randomstr 10`

source /krummellab/data1/ipi/software/cellranger/usr/SOURCE_THIS

mkdir -p /scratch/${USER}/cellranger_mkfastq_${FLOWCELLID}_${uuid}/working  /scratch/${USER}/cellranger_mkfastq_${FLOWCELLID}_${uuid}/working
cd /scratch/${USER}/cellranger_mkfastq_${FLOWCELLID}_${uuid}/working 
trap "{ rm -rf /scratch/${USER}/cellranger_mkfastq_${FLOWCELLID}_${uuid} ; }" EXIT

fastq_dir=$(dirname $(pwd))/fastqs

lanes_argstring=" "
if [ ${LANES} != "ALL" ]
then
    lanes_argstring="--lanes=${LANES} "
fi

echo "running command: "
echo "cellranger-${CELLRANGERVERSION} mkfastq --csv=${SAMPLESHEET} \
                   --run=${BCLDIR} \
                   --barcode-mismatches=${BARCODEMISMATCHES} \
                   ${lanes_argstring} \
                   --output-dir=${fastq_dir} \
                   --localcores=${PBS_NUM_PPN} \
                   --localmem=${MEMORY} "

cellranger-${CELLRANGERVERSION} mkfastq --csv=${SAMPLESHEET} \
                   --run=${BCLDIR} \
                   --barcode-mismatches=${BARCODEMISMATCHES} \
                   ${lanes_argstring} \
                   --output-dir=${fastq_dir} \
                   --localcores=${PBS_NUM_PPN} \
                   --localmem=${MEMORY}


find ${fastq_dir}/ -name "*.gz" -exec dirname {} \; | sort -u | while read subfolder
do
  cd ${subfolder}
  echo "Validating entries within ${subfolder}"
  echo "Testing Gzip integrity"
  ls *gz | xargs -n1 /krummellab/data1/${USER}/aarao_scripts/bash/on_path/test_pigz_integrity
  echo "Generating md5sums"
  md5sum *.fastq.gz > md5checksum.txt
  echo ""
done

echo "Filesizes for All generated fastqs:"
find ${fastq_dir} -name "*.gz" -exec du -h {} \; | sort -k2d


if [ -f ${OUTDIR} ]
then
    mv ${fastq_dir} ${OUTDIR}/${FLOWCELLID}_mkfastq_${uuid}
else
    mkdir -p $(dirname ${OUTDIR})
    mv ${fastq_dir} ${OUTDIR}
fi
