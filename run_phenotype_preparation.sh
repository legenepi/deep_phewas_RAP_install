#!/bin/bash

. RAP.config

INPUTS=${PHENOTYPES_FILTERED}.dx.json

if ! APPLET=`dx ls --brief ${PROJECT_DIR}/phenotype_preparation`; then
    echo "Applet phenotype_preparation not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
elif [ ! -s $INPUTS ]; then
    echo "$INPUTS not found, have you run make_inputs_phenotype_preparation.sh?"
    exit 2
fi

dx mkdir -p $PHENOTYPE_TABLES &&
dx run --destination $PHENOTYPE_TABLES --brief -y -f $INPUTS $APPLET
