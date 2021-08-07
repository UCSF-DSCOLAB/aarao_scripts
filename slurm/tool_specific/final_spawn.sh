#!/bin/bash

set -e
set -o nounset

maxpad=$(echo -e "${GLOBAL_OPTIONAL_HELPTEXT}${REQUIRED_HELPTEXT}${LOCAL_OPTIONAL_HELPTEXT-''}" | grep -v "^#" | cut -f 1 -d ":" | wc -L)

padding=$(printf '%*s' "$maxpad" )

echo "Running using the following options:"

for ht_type in REQUIRED_HELPTEXT LOCAL_OPTIONAL_HELPTEXT GLOBAL_OPTIONAL_HELPTEXT
do
    echo "${ht_type/HELPTEXT/ARGUMENTS}: "
    var_names=($(echo -e "${!ht_type-''}" | grep -v "^#" | cut -f 1 -d " " | xargs))
    if [[ ${var_names-"EMPTY"} == "EMPTY" ]]
    then
        echo ""
        continue
    fi
    for var_name in ${var_names[@]}
    do
        printf "\t%s%s: %s\n" ${var_name} "${padding:${#var_name}}" ${!var_name}
    done
    echo ""
done


GLOBAL_EXPORT_VARS=$(echo -e "${GLOBAL_OPTIONAL_HELPTEXT}" |
                            grep -v "^#" |
                            cut -f 1 -d " " |
                            xargs -n1 |
                            while read i
                            do
                                echo -e "${i}=${!i},\c"
                            done)

if [[ ${LOCAL_EXPORT_VARS-"EMPTY"} == "EMPTY" ]]
then
    LOCAL_EXPORT_VARS=""
fi

sbatch --export="ALL,${GLOBAL_EXPORT_VARS}${LOCAL_EXPORT_VARS}" \
       --error=${LOGDIR}/${JOBNAME}_$(date "+%Y_%m_%d_%H_%M_%S").err \
       --output=${LOGDIR}/${JOBNAME}_$(date "+%Y_%m_%d_%H_%M_%S").out \
       --job-name=${JOBNAME} \
       --partition=${PARTITION} \
       --time=${TIME} \
       --ntasks=${NTASKS} \
       --cpus-per-task ${CPUSPERTASK} \
       --mem-per-cpu ${MEMPERCPU} \
       $(dirname ${0})/run.sh
