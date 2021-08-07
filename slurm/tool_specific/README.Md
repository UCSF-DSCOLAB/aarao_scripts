# Description

This folder provides an easy way to run (mostly) bioinformatic tools of interest on a slurm-powered
cluster using singularity containers. This folder assumes you have closed the entire `aarao_scripts`
repository at the location `/krummellab/data1/${USER}/aarao_scripts`.

# How this directory works

Each subdirectory contains a `spawn.sh` and a `run.sh` file. `spawn.sh` accepts, parses and exports important variables to an sbatch job that calls
the script `run.sh`.

In an attempt to make generalizable spawn/run pairs, every spawn.sh file 
will define, parse and export args required by run.sh but final_spawn.sh will
actually spawn the job. Every spawn.sh script have roughly the following structure.
```
set -e
set -onounset

# Parse the arguments and define global args shared by all
# programs
source $(dirname ${0})/../template_argparse.sh
# Source some essentials if required
source $(dirname ${0})/../../../bash/essentials.sh

##################################
##### This is mandatory
##################################
read -r -d '' REQUIRED_HELPTEXT  << EOF || true
## REQUIRED PARAMETERS ##
REQUIREDARGUMENT1         : Description
...
REQUIREDARGUMENTn         : Description
EOF

##################################
##### This may or may not be there
##################################
read -r -d '' LOCAL_OPTIONAL_HELPTEXT  << EOF || true
## OPTIONAL PARAMETERS ##
OPTIONALARGUMENT1         : Description
...
OPTIONALARGUMENT1         : Description
EOF

##################################
##### This may or may not be there
##################################
read -r -d '' GLOBAL_OVERRIDE_HELPTEXT  << EOF || true
## GLOBAL OVERRIDES ##
GLOBALOVERRIDE1 : Description. 
...
GLOBALOVERRIDEn : Description
EOF

##################################
##### Validate all REQUIREDD and 
##### OPTIONAL args here
##################################

##################################
##### Export local vars so they're 
##### visibe to run.sh
##################################
LOCAL_EXPORT_VARS="\
CONTAINER=$(readlink -e ${CONTAINER}),\
SAMPLE_YML=$(readlink -e ${SAMPLE_YML}),\
WORKING_FOLDER=${WORKING_FOLDER},\
RSCRIPTS_DIR=$(readlink -e ${RSCRIPTS_DIR}),\
GENESET_DIR=$(readlink -e ${GENESET_DIR})"

# A job name for the logs
JOBNAME=process_10x_seurat_$(randomstr 10)

# Call the final spawn script
source $(dirname ${0})/../final_spawn.sh
```
