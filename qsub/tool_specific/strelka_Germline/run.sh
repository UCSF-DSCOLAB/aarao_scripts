#!/bin/bash
set -e
set -o nounset

source /home/${USER}/miniconda3/etc/profile.d/conda.sh
conda activate strelka

mkdir /scratch/${USER}/strekla_GL_${SAMPLE} && cd /scratch/${USER}/strekla_GL_${SAMPLE} 
trap "{ rm -rf /scratch/${USER}/strekla_GL_${SAMPLE} /scratch/${USER}/${SAMPLE}_javatmp ; }" EXIT

extension=${SAMFILE##*.}

out_dir=$(dirname ${SAMFILE})

/krummellab/data1/ipi/software/strelka/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${SAMFILE} \
    --referenceFasta ${GENOMEREF} \
    --exome \
    --runDir ${SAMPLE}

# This removes the stupid smtp issue that hangs the tool
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${SAMPLE}/runWorkflow.py 

${SAMPLE}/runWorkflow.py \
       --mode local \
       --jobs ${PBS_NUM_PPN} \
       --memGb ${MEMORY}

mv /scratch/${USER}/strekla_GL_${SAMPLE}/${SAMPLE}/results/ ${out_dir}/strelka/
