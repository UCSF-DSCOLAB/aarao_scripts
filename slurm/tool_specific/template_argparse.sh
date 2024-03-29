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
LOGDIR             : A path to a logdir (default: /krummellab/data1/${USER}/logs)
PARTITION          : Slurm partition to use (default: krummellab)
TIME               : Max job runtime (default: 1-00:00:00)
NTASKS             : Number of tasks in this job (default: 1) [YMMV with modifying this parameter!]
CPUSPERTASK        : Cpus required per task (default: 12)
MEMPERCPU          : Memory required per cpu (default: 10gb)
MEMORYFRAC         : Fraction of requested memory usable by the task (default: 0.9) 
USABLEMEMORY       : Memory usable by the task (default: CPUSPERTASK x MEMPERCPU x MEMORYFRAC)
COLLAPSEDIRSCRIPT  : /krummellab/data1/${USER}/aarao_scripts/python/collapse_folderlist.py
ECHO_CMD           : Positional argument to echo the sbatch command (Useful for debugging exported vars)
DRY_RUN            : Positional argument (Useful for debugging exported vars, implies ECHO_CMD)
EOF

DEFAULTLOGDIR=/krummellab/data1/${USER}/logs
DEFAULTPARTITION=krummellab
DEFAULTTIME=1-00:00:00
DEFAULTNTASKS=1
DEFAULTCPUSPERTASK=12
DEFAULTMEMPERCPU=10gb
DEFAULTMEMORYFRAC=0.9
DEFAULTCOLLAPSEDIRSCRIPT=/krummellab/data1/${USER}/aarao_scripts/python/collapse_folderlist.py

GLOBAL_SKIP_DEFAULT_ARGS=("USABLEMEMORY" "ECHO_CMD" "DRY_RUN")

# POSITIONAL_ARGS can be empty and set -o nounset will raies an error
if [[ ${POSITIONAL_ARGS[@]-"EMPTY"} == "EMPTY" ]]
then
    POSITIONAL_ARGS=""
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
  if [[ ${GLOBAL_OVERRIDE_HELPTEXT-"EMPTY"} != "EMPTY" ]]
  then
    echo -e "${GLOBAL_OVERRIDE_HELPTEXT-''}\n"
  fi
  echo -e "${GLOBAL_OPTIONAL_HELPTEXT}\n"
  exit 1
}
