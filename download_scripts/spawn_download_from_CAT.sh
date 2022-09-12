#!/bin/bash
set -e

ALLOWED_RUN_MODES=(qsub slurm direct)
if [ "$(pwd)" != "/krummellab/data1/immunox/staging" ]
then
  echo "Please only run this script from /krummellab/data1/immunox/staging"
  exit
fi
if [ $# -ne 2 ]
then
  echo "Need a run_mode, and one file to download data from"
  echo "USAGE: spawn_download_from_CAT.sh <run_mode> <download_manifest>"
  exit
else
    if ! [[ ${ALLOWED_RUN_MODES[@]} =~ $1 ]]
    then
        echo "run_mode must be one of the following: ${ALLOWED_RUN_MODES[@]}"
        exit
    fi
    filename=`readlink -e $2 || echo "NOFILE"`
    if [[ "${filename}" == "NOFILE" ]]
    then
        echo "provided filename doesn't exist"
        exit
    fi
fi

while read line
do
  if [[ $line != */ ]]
  then
    echo "ERROR: A line in the input file did not end with a /"
    exit 1
  fi
done < ${filename}

if [[ ! -f ~/.netrc ]]
then
  echo "ERROR: Cannot find a credentials file. This script reads credentials from ~/.netrc."
  exit 1
elif ! grep -q "fastq.ucsf.edu" ~/.netrc
then
  echo "ERROR: ~/netrc does not contain a record for fastq.ucsf.edu. Cannot continue without proper credentials for that server."
  echo "A proper netrc record looks like:"
  echo """
machine fastq.ucsf.edu
        login hiseq_user
        password <PASSWORD>
"""
  exit 1
fi


identifier=$(basename ${filename%.*})
log_out=/krummellab/data1/${USER}/logs/${identifier}_$(date "+%Y_%m_%d_%H_%M_%S").out
log_err=/krummellab/data1/${USER}/logs/${identifier}_$(date "+%Y_%m_%d_%H_%M_%S").err
if [[ "${1}" == "direct" ]]
then
    export FILENAME=${filename}
    export TMPDIR=/scratch/${USER}
    bash download_from_CAT.sh 1>${log_out} 2>${log_err}
elif [[ "${1}" == "slurm" ]]
then
    echo NotImplemented
elif [[ "${1}" == "qsub" ]]
then
    qsub -v "FILENAME=${filename}" \
         -e ${log_err} \
         -o ${log_out} \
         -N ${identifier} \
         -l nodes=1:ppn=1 \
         -l vmem=10G \
         -d $(dirname $0) \
         -w $(dirname $0) \
         download_from_CAT.sh
else
    echo "ERROR"
fi
