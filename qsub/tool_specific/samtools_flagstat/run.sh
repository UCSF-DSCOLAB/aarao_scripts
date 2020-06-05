#!/bin/bash
set -e
set -o nounset

source /data/shared/krummellab/ipi/ipi_usr/SOURCE_THIS 

samtools flagstat ${BAMFILE} > ${BAMFILE%.bam}.flagstat