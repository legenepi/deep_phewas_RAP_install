#!/bin/bash

. RAP.config

for i in $phenotype_files $kinship_file; do
    if ! dx ls $i > /dev/null; then
        echo "Required intermediate file $i missing from ${PHENOTYPES}, have you run phenotype_generation.sh?"
        exit 1
    fi
done

FILTER_FILES="groupings phewas_manifest"
STRING_VARS="phenotype_filtered_save_name stats_save"
BOOL_VARS="relate_remove IVNT"

export $FILTER_FILES kinship_file $STRING_VARS $BOOL_VARS

Rscript - <<-RSCRIPT
    suppressMessages(library(tidyverse))
    suppressMessages(library(jsonlite))
    source("R/make_inputs_functions.R")
    
    phenotype_files <- list(phenotype_files="$phenotype_files" %>%
        str_split_1(" ") %>%
        map(get_file_id))
    kinship_file <- get_config("kinship_file") %>%
        map(get_file_id)
    filter_files <- get_config("$FILTER_FILES") %>%
        map(~get_upload_id(., "$PROJECT_ID", "$PROJECT_DIR"))
    string_vars <- get_config("$STRING_VARS")
    bool_vars <- get_config("$BOOL_VARS") %>%
        map(as.logical)
    c(phenotype_files, kinship_file, filter_files, string_vars, bool_vars) %>%
        write_json("${PHENOTYPES_FILTERED}.json", pretty=TRUE, auto_unbox=TRUE)
RSCRIPT

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER
java -jar $DXCOMPILER compile WDL/phenotype_preparation.wdl -project $PROJECT_ID -compileMode IR -inputs ${PHENOTYPES_FILTERED}.json
