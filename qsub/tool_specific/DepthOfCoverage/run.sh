#!/bin/bash
set -e
set -o nounset
module load CBC gatk/${GATKVERSION}

mkdir /scratch/${USER}/DepthOfCoverage_${SAMPLE}
cd /scratch/${USER}/DepthOfCoverage_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/DepthOfCoverage_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT

out_dir=$(dirname ${SAMFILE})

if [ ${INTERVALFILE-"EMPTY"} == "NONE" ]
then
    INTERVALSTRING=" "
else
    INTERVALSTRING="--intervals ${INTERVALFILE}"
fi

java -Djava.io.tmpdir=/scratch/${USER}/${SAMPLE}_javatmp \
    -Xmx${MEMORY}g \
    -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -I ${SAMFILE} \
    ${INTERVALSTRING} \
    --reference_sequence ${GENOMEREF} \
    --out ${SAMPLE} \
    -nt ${PBS_NUM_PPN} \
    --omitIntervalStatistics \
    --omitDepthOutputAtEachBase \
    --summaryCoverageThreshold 10 \
    --summaryCoverageThreshold 25 \
    --summaryCoverageThreshold 50 \
    --summaryCoverageThreshold 100 

mv /scratch/${USER}/DepthOfCoverage_${SAMPLE}/ ${out_dir}/DepthOfCoverage/
