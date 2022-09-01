#!/bin/bash

if [[ ! ${HOSTNAME} =~ ^c4 ]]
then
    source /krummellab/data1/arrao/software/lftp/SOURCE_THIS
fi

set -e
set -o nounset

prefix=$(basename ${FILENAME%.*})
mkdir ${TMPDIR}/${prefix}/
cd ${TMPDIR}/${prefix}/
trap "{ rm -rf ${TMPDIR}/${prefix}/ ; }" EXIT


while read folder
do
    echo "mirror -P 10 ${folder}" >> lftp_commands.txt
done < ${FILENAME}

lftp -p 22 --user hiseq_user sftp://fastq.ucsf.edu < lftp_commands.txt
rsync -a ${TMPDIR}/${prefix}/ /krummellab/data1/immunox/staging/from_CAT/${prefix}_$(date "+%Y_%m_%d_%H_%M_%S")/
