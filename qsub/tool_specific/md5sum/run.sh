#!/bin/bash

set -e
set -o nounset

cd ${INPUT_DIR} && \
    /usr/bin/md5sum * > md5checksum.txt
