#!/bin/bash

set -e
set -o nounset

for gov in $(echo -e "${GLOBAL_OPTIONAL_HELPTEXT}" | grep -v "^#" | cut -f 1 -d " " | xargs)
do
    if [[ "${GLOBAL_SKIP_DEFAULT_ARGS[@]}" =~ "${gov}" ]]
    then
        continue
    elif [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "${gov}" ]]
    then
        temp=DEFAULT${gov}
        declare "${gov}=${!temp}"
    fi
done

if [[ ! "${RECEIVED_NAMED_ARGS[@]}" =~ "USABLEMEMORY" ]]
then
    USABLEMEMORY=$(echo "${MEMPERCPU} * ${CPUSPERTASK} * ${MEMORYFRAC} / 1"  | sed s/gb//g | bc)
fi

if [[ "${POSITIONAL_ARGS[@]}" =~ "ECHO_CMD" ]]
then
    ECHO_CMD=TRUE
else
    ECHO_CMD=FALSE
fi


if [[ "${POSITIONAL_ARGS[@]}" =~ "DRY_RUN" ]]
then
    ECHO_CMD=TRUE
    DRY_RUN=TRUE
else
    DRY_RUN=FALSE
fi

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


if [[ "${ECHO_CMD}" == "TRUE" ]]
then
    echo sbatch --export="ALL,${GLOBAL_EXPORT_VARS}${LOCAL_EXPORT_VARS}" \
                --error=${LOGDIR}/${JOBNAME}_$(date "+%Y_%m_%d_%H_%M_%S").err \
                --output=${LOGDIR}/${JOBNAME}_$(date "+%Y_%m_%d_%H_%M_%S").out \
                --job-name=${JOBNAME} \
                --partition=${PARTITION} \
                --time=${TIME} \
                --ntasks=${NTASKS} \
                --cpus-per-task ${CPUSPERTASK} \
                --mem-per-cpu ${MEMPERCPU} \
                $(dirname ${0})/run.sh
fi

if [[ "${DRY_RUN}" == "TRUE" ]]
then
    echo -e "\n\`DRY_RUN\` was requested.... exiting."
    exit 1
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
