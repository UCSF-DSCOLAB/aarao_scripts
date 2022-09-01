#!/bin/bash

set -e
set -o nounset

prefix=$(basename ${FILENAME%.*})
mkdir ${TMPDIR}/${prefix}/
cd ${TMPDIR}/${prefix}/
trap "{ rm -rf ${TMPDIR}/${prefix}/ ; }" EXIT

export WGETRC=/krummellab/data1/immunox/staging/wgetrcs/${IHGUSERNAME}
wget --no-check-certificate \
     -m \
     -R "index.html*" \
     --no-parent \
     -c \
     -i ${FILENAME}

if [ $? -ne 0 ]
then
    echo "FAILED"
    exit 1
fi

cwd=$(pwd)
while read entry
do
  cd ${cwd}${entry#https:/}
  find $(pwd) -name "*.gz" -exec dirname {} \; | sort -u | while read subfolder
  do
    cd ${subfolder}
    echo "Validating entries within ${subfolder}"
    echo "#### Gzip integrity ####"
    ls *gz | xargs -n1 /krummellab/data1/${USER}/scripts/bash/on_path/test_gzip_integrity
    echo "#### md5sums ####"
    if [ -f md5checksum.txt ]
    then
        echo "Validating md5sums for entry ${entry}"
        md5sum --check md5checksum.txt
    else
       echo "No md5checksum.txt file for this entry"
    fi
  done
done < ${FILENAME}

rsync -a ${TMPDIR}/${prefix}/ihg-client.ucsf.edu/ /krummellab/data1/immunox/staging/ihg-client.ucsf.edu/

