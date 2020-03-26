#!/bin/bash
set -e
set -o nounset

source /krummellab/data1/ipi/ipi_usr/SOURCE_THIS

mkdir /scratch/arrao/fastp_${SAMPLE} && cd /scratch/arrao/fastp_${SAMPLE} 
trap "{ rm -rf /scratch/arrao/fastp_${SAMPLE} ; }" EXIT

getExtension() {
      x=${1##*.}
      if [ ${x} == "gz" ]
      then
            echo .`getExtension ${1%.*}`.${x}
      elif [ ${x} == "fastq" ]
      then
            echo ${x}
      elif [ ${x} == "fastq" ]
      then
            echo ${x}
      else
            echo "Could not identify extension"
            exit 1
      fi
}

extension=`getExtension ${FQ1}`

FQ2argstring=" "
if [ $(basename ${FQ2}) != "EMPTY" ]
then
      FQ2argstring="-I ${FQ2} \
                    -O `basename ${FQ2%${extension}}`_trimmed${extension} \
                    --adapter_sequence_r2 CTGTCTCTTATACACATCT \
                    --correction "
fi

fastp -i ${FQ1} \
      -o `basename ${FQ1%${extension}}`_trimmed${extension} \
      --adapter_sequence CTGTCTCTTATACACATCT \
      ${FQ2argstring} \
      --length_required 20 \
      --trim_poly_g  \
      --thread ${PBS_NUM_PPN} \
      -j ${SAMPLE}_fastp.json \
      -h ${SAMPLE}_fastp.html

mkdir -p ${OUTDIR}
mv /scratch/arrao/fastp_${SAMPLE}/* ${OUTDIR}/
