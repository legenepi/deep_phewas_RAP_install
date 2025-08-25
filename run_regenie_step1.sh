#!/bin/bash

. RAP.config

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/regenie_step1`; then
    echo "Workflow assocition_testing not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
fi

if [ ! -s ${REGENIE_STEP1_INPUTS}.dx.json ]; then
    echo "Regenie step1 inputs specification file ${REGENIE_STEP1_INPUTS}.dx.json not found, have you run make_inputs_regenie_step1.sh?"
    exit 2
fi

dx run --ssh --debug-on All --destination ${PROJECT_DIR}/step1 --brief -y -f  ${REGENIE_STEP1_INPUTS}.dx.json $WORKFLOW
