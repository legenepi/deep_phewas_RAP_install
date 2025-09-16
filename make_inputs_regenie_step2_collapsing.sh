#!/bin/bash

. RAP.config

if ! WORKFLOW=`dx ls --brief ${PROJECT_DIR}/regenie_step2_collapsing`; then
    echo "Workflow regenie_step2_collapsing not found in ${PROJECT_DIR} on RAP, have you run install_workflows.sh?"
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


# "regenie_step2_collapsing.mask": "File",
# "regenie_step2_collapsing.loco_bt": "Array[Array[File]]",
# "regenie_step2_collapsing.regenie_bt.extract": "File? (optional)",
# "regenie_step2_collapsing.phenotypes_bt": "Array[File]",
# "regenie_step2_collapsing.pred_list_qt": "Array[File]",
# "regenie_step2_collapsing.covariates": "File? (optional)",
# "regenie_step2_collapsing.aaf_bins": "Array[Float]",
# "regenie_step2_collapsing.catCovarList": "String? (optional)",
# "regenie_step2_collapsing.analysis_name": "String",
# "regenie_step2_collapsing.setlist": "File",
# "regenie_step2_collapsing.tests": "String",
# "regenie_step2_collapsing.regenie_qt.extract": "File? (optional)",
# "regenie_step2_collapsing.maxaff": "Float",
# "regenie_step2_collapsing.phenotypes_qt": "Array[File]",
# "regenie_step2_collapsing.pred_list_bt": "Array[File]",
# "regenie_step2_collapsing.genos": "Array[Array[File]]",
# "regenie_step2_collapsing.loco_qt": "Array[Array[File]]",
# "regenie_step2_collapsing.annot": "File",
# "regenie_step2_collapsing.phenoColList": "String? (optional)",
# "regenie_step2_collapsing.exclude": "File? (optional)",
# "regenie_step2_collapsing.covarColList": "String? (optional)",
# "regenie_step2_collapsing.joint_tests": "String",
# "regenie_step2_collapsing.gene_list": "File"

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

    chroms <- str_split_1("$SET_CHROMS", " ")

    snp_list_rsid <- "$snp_list" %>%
        str_replace(".csv", "_rsid.csv")

    snp_list %>%
        pull(rsid) %>%
        cat(file=snp_list_rsid, sep="\n")

    list(regenie_step2_collapsing.phenotypes_bt = phenotypes_bt %>% map(get_file_id),
            regenie_step2_collapsing.phenotypes_qt = phenotypes_qt %>% map(get_file_id),
            regenie_step2_collapsing.genos = get_genos("$EXOME_PATH", chroms=chroms, format="plink"),
            regenie_step2_collapsing.covariates = get_upload_id("$covar", "$PROJECT_DIR"),
            regenie_step2_collapsing.covarColList = "$covarColList",
            regenie_step2_collapsing.catCovarList = "$catCovarList",
            regenie_step2_collapsing.phenoColList = "$phenoColList",
            regenie_step2_collapsing.snp_list = get_upload_id(snp_list_rsid, "$PROJECT_DIR"),
            regenie_step2_collapsing.pred_list_bt = pred_files_bt %>% map(get_file_id),
            regenie_step2_collapsing.loco_bt = pred_files_bt %>% map(~get_loco(., "${PROJECT_DIR}/step1")),
            regenie_step2_collapsing.pred_list_qt = pred_files_qt %>% map(get_file_id),
            regenie_step2_collapsing.loco_qt = pred_files_qt %>% map(~get_loco(., "${PROJECT_DIR}/step1")),
            regenie_step2_collapsing.annot = get_file_id("$ANNOT"),
            regenie_step2_collapsing.setlist = get_file_id("$SETS"),
            regenie_step2_collapsing.tests = "$REGENIE_STEP2_TESTS",
            regenie_step2_collapsing.maxaff = $REGENIE_STEP2_MAXAFF,
            regenie_step2_collapsing.joint_tests = "$REGENIE_STEP2_JOINT_TESTS",
            regenie_step2_collapsing.aaf_bins = "$REGENIE_STEP2_AAF_BINS",
            regenie_step2_collapsing.analysis_name = "$analysis_name",
            regenie_step2_collapsing.gene_list = get_upload_id("$gene_list", "$PROJECT_DIR") %>%
        write_json("${REGENIE_STEP2_INPUTS}.json", pretty=TRUE, auto_unbox=TRUE)
RSCRIPT

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER

java -jar $DXCOMPILER compile WDL/regenie_step2_collapsing.wdl -project $PROJECT_ID -compileMode IR -inputs ${REGENIE_STEP2_INPUTS}.json
