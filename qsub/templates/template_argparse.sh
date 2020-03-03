#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../../../bash/argparse.sh


USAGE_HELPTEXT="USAGE:  spawn.sh REQD_ARG1=RVAL1 ... REQD_ARGn=RVALn OPTNL_ARG1=OVAL1 ... OPTNL_ARGn=OVALn POSITIONAL_ARGS"

read -r -d '' GLOBAL_OPTIONAL_HELPTEXT << EOF || true
## GLOBAL OPTIONAL PARAMETERS ##
LOGDIR=/path/to/logdir (/krummellab/data1/arrao/logs)
NODEREQS=nodes=1:ppn=64
MEMREQS=vmem=500gb
EOF

if [ ${LOGDIR-"EMPTY"} == "EMPTY" ]
then
    LOGDIR=/krummellab/data1/arrao/logs
else
    LOGDIR=${LOGDIR#/}
fi

if [ ${NODEREQS-"EMPTY"} == "EMPTY" ]
then
    NODEREQS=nodes=1:ppn=64
fi

if [ ${MEMREQS-"EMPTY"} == "EMPTY" ]
then
    MEMREQS=vmem=500gb
fi


function print_help() {
  echo -e "${USAGE_HELPTEXT}"
  echo -e "${REQUIRED_HELPTEXT}"
  if [ ${LOCAL_OPTIONAL_HELPTEXT-"EMPTY"} != "EMPTY" ]
  then
    echo -e "${LOCAL_OPTIONAL_HELPTEXT-''}"
  fi
  echo -e "${GLOBAL_OPTIONAL_HELPTEXT}"
  exit 1
}
