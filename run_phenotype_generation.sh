#!/bin/bash

. RAP.config

INPUTS=${PHENOTYPES}.dx.json

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/phenotype_generation`; then
    echo "Workflow phenotype_generation not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
elif [ ! -s $INPUTS ]; then
    echo "$INPUTS not found, have you run make_inputs_phenotype_generation.wdl?"
    exit 2
fi

dx mkdir -p $INTERMEDIATE &&
dx run --destination $INTERMEDIATE --brief -y -f $INPUTS $WORKFLOW
