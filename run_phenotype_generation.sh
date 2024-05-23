#!/bin/bash

. RAP.config

INPUTS=${PHENOTYPES_GENERATED}.dx.json

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/phenotype_generation`; then
    echo "Workflow phenotype_generation not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
elif [ ! -s $INPUTS ]; then
    echo "$INPUTS not found, have you run make_inputs_phenotype_generation.sh?"
    exit 2
fi

dx mkdir -p $PHENOTYPES &&
dx run --destination $PHENOTYPES --brief -y -f $INPUTS $WORKFLOW
