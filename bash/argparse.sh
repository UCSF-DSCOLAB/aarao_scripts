#!/bin/bash
set -e
set -o nounset

POSITIONAL_ARGS=()
RECEIVED_NAMED_ARGS=()
for i in "$@"
do
case $i in
    *=*)
    declare "${i%%=*}"="${i#*=}"
    RECEIVED_NAMED_ARGS+=("${i%=*}")
    shift # past argument=value
    ;;
    *)
    POSITIONAL_ARGS+=("${i}")
    ;;
esac
done

echo "ARGPARGE: RECEIVED_NAMED_ARGS="${RECEIVED_NAMED_ARGS[@]-""}
echo "ARGPARSE: POSITIONAL_ARGS="${POSITIONAL_ARGS[@]-""}
echo -e "\n"
