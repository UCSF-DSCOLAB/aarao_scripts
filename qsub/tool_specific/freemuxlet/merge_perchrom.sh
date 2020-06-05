#!/bin/bash
source /krummellab/data1/ipi/software/popscle/usr/SOURCE_THIS

set -e
set -o nounset

if [ ! -f "${OUTDIR}/expected_chroms.list" ]
then
    echo "Cannot find expected_chroms.list"
fi

prefix=$(basename ${BAMFILE%.*})
extension=${BAMFILE##*.}

timer=0
while :
do
    sleep ${SLEEPTIME}
    files_remaining=0
    timer=$((timer+SLEEPTIME))
    if [ ${timer} -gt ${TIMEOUT} ]
    then
        echo 'Timed out waiting for the tasks to finish!'
        exit
    fi

    while read chrom
    do
        if [ ! -f "${OUTDIR}/${prefix}_${chrom}.cel.gz" ]
        then
            files_remaining=$((files_remaining+1))
        fi
    done < ${OUTDIR}/expected_chroms.list

    if [ ${files_remaining} -eq 0 ]
    then
        echo 'All files received. Let us continue!'
        break
    else
        echo "Waiting for ${files_remaining} more files"
    fi
done

~/miniconda3/envs/pysam/bin/python /krummellab/data1/arrao/scripts/python/merge_perchrom_dsc_pileups.py \
    ${OUTDIR}/${prefix}_ ${OUTDIR}/${prefix}_MERGED ${BAMFILE} 

popscle freemuxlet --plp ${OUTDIR}/${prefix}_MERGED \
                   --out ${OUTDIR}/${prefix} \
                   --nsample ${NUMSAMPLES} \
                   --seed ${RANDOMSEED} \
                   --group-list ${BARCODELIST}
