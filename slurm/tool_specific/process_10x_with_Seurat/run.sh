#!/bin/bash

set -e
set -o nounset

if [[ ${TMPDIR-"EMPTY"} ==  "EMPTY" ]]
then
    TMPDIR=/scratch/${USER}/$(head -10 /dev/null | md5sum | cut -c 1-10)
    trap "{ rm -rf ${TMPDIR} ; }" EXIT
fi

if [ "${WORKING_FOLDER}" == "TMPDIR" ]
then
  WORKING_FOLDER=${TMPDIR}
fi

yml_folders=($(grep "datadir\|OUT_FOLDER\|freemuxlet" ${SAMPLE_YML}  | 
                        xargs -n2 |
                        cut -f 2 -d " " |
                        xargs -n1 dirname | 
                        sort -u
                        ))

bindmount_string=$(python3 ${COLLAPSEDIRSCRIPT} --prefixB $(dirname ${SAMPLE_YML}) ${yml_folders[@]} ${RSCRIPTS_DIR} ${GENESET_DIR} ${WORKING_FOLDER})

singularity exec \
      ${bindmount_string[@]} \
      -B ${TMPDIR}:/tmp/ \
      --pwd ${WORKING_FOLDER} \
      ${CONTAINER} Rscript \
            ${RSCRIPTS_DIR}/process_10x_with_Seurat.R \
                SAMPLE_YML=${SAMPLE_YML} \
                RSCRIPTS_DIR=${RSCRIPTS_DIR} \
                GENESET_DIR=${GENESET_DIR} \
                WORKING_FOLDER=${WORKING_FOLDER}
