#!/bin/bash

. RAP.config

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/regenie_step1`; then
    echo "Workflow assocition_testing not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
fi

if ! PHENOTYPES=`dx ls --brief "${PHENOTYPE_TABLES}/*${phenotype_filtered_save_name}*"`; then
    echo "Phenotype files not found, have you run phenotype_preparation.sh"
    exit 2
fi

Rscript - <<-RSCRIPT
    suppressMessages(library(tidyverse))
    suppressMessages(library(jsonlite))
    source("R/make_inputs_functions.R")

    groupings <- read_tsv("$groupings") %>%
        pull(group) %>%
        unique

    phenotypes <- paste0("${PROJECT_ID}:${PHENOTYPE_TABLES}/", groupings, "_", "${phenotype_filtered_save_name}.gz") 

    list(regenie_step1.genos = get_genos("$STEP1_GENO", chroms=as.character(1:22)),
            regenie_step1.missing_thresh = $REGENIE_STEP1_MISSING_THRESH,
            regenie_step1.prepared_phenotypes = phenotypes %>% map(get_file_id),
            regenie_step1.covar = get_upload_id("$covar", "$PROJECT_DIR"),
            regenie_step1.covarColList = "$covarColList",
            regenie_step1.catCovarList = "$catCovarList",
            regenie_step1.phewas_manifest = get_upload_id("$phewas_manifest", "$PROJECT_DIR"),
            regenie_step1.phenoColList = "$phenoColList") %>%
        write_json("${REGENIE_STEP1_INPUTS}.json", pretty=TRUE, auto_unbox=TRUE)
RSCRIPT

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER

java -jar $DXCOMPILER compile WDL/regenie_step1.wdl -project $PROJECT_ID -compileMode IR -inputs ${REGENIE_STEP1_INPUTS}.json
