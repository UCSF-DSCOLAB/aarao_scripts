#!/bin/bash
set -e
set -o nounset

source $(dirname ${0})/../../../bash/argparse.sh


read -r -d '' USAGE_HELPTEXT << EOF || true
USAGE:  spawn.sh REQD_ARG1=RVAL1 \\
                 REQD_ARG2=RVAL2 \\
                 REQD_ARGn=RVALn \\
                 OPTNL_ARG1=OVAL1 \\
                 OPTNL_ARG2=OVAL2 \\
                 OPTNL_ARGn=OVALn \\
                 POSITIONAL_ARG1 \\
                 POSITIONAL_ARG2 \\
                 POSITIONAL_ARGn
EOF


read -r -d '' GLOBAL_OPTIONAL_HELPTEXT << EOF || true
## GLOBAL OPTIONAL PARAMETERS ##
LOGDIR=/path/to/logdir (/krummellab/data1/${USER}/logs)
NODEREQS=nodes=1:ppn=64
MEMREQS=vmem=500gb
EOF

if [ ${LOGDIR-"EMPTY"} == "EMPTY" ]
then
    LOGDIR=/krummellab/data1/${USER}/logs
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
  echo "################################################################################"
  echo "################################################################################"
  echo -e "${USAGE_HELPTEXT}\n"
  echo -e "${REQUIRED_HELPTEXT}\n"
  if [[ ${LOCAL_OPTIONAL_HELPTEXT-"EMPTY"} != "EMPTY" ]]
  then
    echo -e "${LOCAL_OPTIONAL_HELPTEXT-''}\n"
  fi
  echo -e "${GLOBAL_OPTIONAL_HELPTEXT}\n"
  exit 1
}
