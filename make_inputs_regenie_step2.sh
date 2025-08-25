#!/bin/bash

. RAP.config

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/regenie_step2_SV`; then
    echo "Workflow assocition_testing not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
    exit 1
fi

readarray pheno_files_bt < <(awk 'NR > 1 { group[$2] } END { for (i in group) print i"_'${phenotype_filtered_save_name}'_bin.txt" }' $groupings)
readarray pheno_files_qt < <(awk 'NR > 1 { group[$2] } END { for (i in group) print i"_'${phenotype_filtered_save_name}'_quant.txt" }' $groupings)

readarray pred_files_bt < <(awk 'NR > 1 { group[$2] } END { for (i in group) print i"_'${phenotype_filtered_save_name}'_bin_step1_pred.list" }' $groupings)
readarray pred_files_qt < <(awk 'NR > 1 { group[$2] } END { for (i in group) print i"_'${phenotype_filtered_save_name}'_quant_step1_pred.list" }' $groupings)

for i in ${pheno_files_bt[@]} ${pheno_files_qt[@]} ${pred_files_bt[@]} ${pred_files_qt[@]}; do
    if ! dx ls ${PROJECT_DIR}/step1/$i > /dev/null; then
        echo "Step 1 file $i missing, have you run regenie step 1?"
        exit 2
    fi
done


# "regenie_step2_SV.phenotypes_qt": "Array[File]",
# "regenie_step2_SV.loco_qt": "Array[Array[File]]",
# "regenie_step2_SV.loco_bt": "Array[Array[File]]",
# "regenie_step2_SV.phenotypes_bt": "Array[File]",
# "regenie_step2_SV.pred_list_qt": "Array[File]",
# "regenie_step2_SV.covariates": "File? (optional)",
# "regenie_step2_SV.snp_list": "File",
# "regenie_step2_SV.catCovarList": "String? (optional)",
# "regenie_step2_SV.analysis_name": "String",
# "regenie_step2_SV.covarColList": "String? (optional)",
# "regenie_step2_SV.pred_list_bt": "Array[File]",
# "regenie_step2_SV.genos": "Array[Array[File]]",
# "regenie_step2_SV.phenoColList": "String? (optional)",
# "regenie_step2_SV.exclude": "File? (optional)"

Rscript - <<-RSCRIPT
    suppressMessages(library(tidyverse))
    suppressMessages(library(jsonlite))
    source("R/make_inputs_functions.R")

    groupings <- read_tsv("$groupings") %>%
        pull(group) %>%
        unique

    split_env_array <- function(x) {
        str_split_1(x, " ") %>%
        str_remove("\n") %>%
        paste0("${PROJECT_DIR}/step1/", .)
    }
       
    phenotypes_bt <- split_env_array("${pheno_files_bt[@]}")
    phenotypes_qt <- split_env_array("${pheno_files_qt[@]}")
    pred_files_bt <- split_env_array("${pred_files_bt[@]}")
    pred_files_qt <- split_env_array("${pred_files_qt[@]}")

    snp_list <- read_csv("$snp_list")
    chroms <- snp_list %>%
        pull(chromosome) %>%
        unique
    
    snp_list_rsid <- "$snp_list" %>%
        str_replace(".csv", "_rsid.csv")

    snp_list %>%
        pull(rsid) %>%
        cat(file=snp_list_rsid, sep="\n")

    list(regenie_step2_SV.phenotypes_bt = phenotypes_bt %>% map(get_file_id),
            regenie_step2_SV.phenotypes_qt = phenotypes_qt %>% map(get_file_id),
            regenie_step2_SV.genos = get_genos("$STEP2_GENO", chroms=chroms, format="bgen"),
            regenie_step2_SV.covariates = get_upload_id("$covar", "$PROJECT_DIR"),
            regenie_step2_SV.covarColList = "$covarColList",
            regenie_step2_SV.catCovarList = "$catCovarList",
            regenie_step2_SV.phenoColList = "$phenoColList",
            regenie_step2_SV.snp_list = get_upload_id(snp_list_rsid, "$PROJECT_DIR"),
            regenie_step2_SV.pred_list_bt = pred_files_bt %>% map(get_file_id),
            regenie_step2_SV.loco_bt = pred_files_bt %>% map(~get_loco(., "${PROJECT_DIR}/step1")),
            regenie_step2_SV.pred_list_qt = pred_files_qt %>% map(get_file_id),
            regenie_step2_SV.loco_qt = pred_files_qt %>% map(~get_loco(., "${PROJECT_DIR}/step1")),
            regenie_step2_SV.analysis_name = "$analysis_name") %>%
        write_json("${REGENIE_STEP2_INPUTS}.json", pretty=TRUE, auto_unbox=TRUE)
RSCRIPT

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER

java -jar $DXCOMPILER compile WDL/regenie_step2_SV.wdl -project $PROJECT_ID -compileMode IR -inputs ${REGENIE_STEP2_INPUTS}.json
