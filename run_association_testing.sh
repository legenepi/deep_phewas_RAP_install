#!/bin/bash

. RAP.config

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/association_testing`; then
    echo "Workflow assocition_testing not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
fi

if ! PHENOTYPES=`dx ls --brief "${PHENOTYPE_TABLES}/*${phenotype_filtered_save_name}*"`; then
    echo "Phenotype files not found, have you run phenotype_preparation.sh"
    exit 2
fi

dx mkdir -p $RESULTS 

#{
#  "association_testing.bgis": "Array[File]+",
#  "association_testing.phenotypes": "Array[File]+",
#  "association_testing.sample_files": "Array[File]+",
#  "association_testing.snp_list": "File",
#  "association_testing.phenotype_inclusion_file": "File? (optional)",
#  "association_testing.analysis_name": "String",
#  "association_testing.phewas_manifest": "File? (optional)",
#  "association_testing.covariates": "File? (optional)",
#  "association_testing.bgens": "Array[File]+"
#}

OPTIONAL_INPUTS="covariates phenotype_inclusion_file phewas_manifest"
export $OPTIONAL_INPUTS

Rscript - <<-RSCRIPT
    suppressMessages(library(tidyverse))
    suppressMessages(library(jsonlite))
    source("R/make_inputs_functions.R")

    snps <- read_csv("$snp_list")

    genos <- snps %>%
        pull(chromosome) %>%
        unique %>%
        as.character() %>%
        map(~paste0("$PROJECT_ID", ":", str_replace("$IMPUTED", "CHROM", .)))
    
    inputs <- list(
        phenotypes="$PHENOTYPES" %>%
        str_split_1("\n") %>%
        map(~paste0("dx://", .)),
        bgens=paste0(genos, ".bgen") %>% map(get_file_id),
        bgis=paste0(genos, ".bgen.bgi") %>% map(get_file_id),
        sample_files=paste0(genos, ".sample") %>% map(get_file_id),
        snp_list=get_upload_id("$snp_list", "$PROJECT_ID", "$RESULTS"),
        analysis_name="$analysis_name") %>%
        set_names(~paste0("association_testing.", .))
    
    optional_inputs <- get_config("$OPTIONAL_INPUTS", "association_testing") %>%
        map(~get_upload_id(., "$PROJECT_ID", "$RESULTS")) 
    
    c(inputs, optional_inputs) %>%
        write_json("${ASSOC_INPUTS}.json", pretty=TRUE, auto_unbox=TRUE)
RSCRIPT

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER

java -jar $DXCOMPILER compile WDL/association_testing.wdl -project $PROJECT_ID -compileMode IR -inputs ${ASSOC_INPUTS}.json &&
dx run --destination $RESULTS --brief -y -f ${ASSOC_INPUTS}.dx.json $WORKFLOW
