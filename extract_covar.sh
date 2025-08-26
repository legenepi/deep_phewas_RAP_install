#!/bin/bash

. RAP.config 

dx ls ${INPUTS}/$COVAR_FIELDS 2> /dev/null && dx rm -a ${INPUTS}/$COVAR_FIELDS
dx upload --path ${INPUTS}/ $COVAR_FIELDS

dx run table-exporter \
    -idataset_or_cohort_or_dashboard=$DATASET \
    -ientity="participant" \
    -ifield_names_file_txt="${INPUTS}/$COVAR_FIELDS" \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ioutput="$COVAR" \
    --destination $INPUTS \
    --instance-type mem2_ssd2_x16 \
    --name extract_covar \
    --ignore-reuse \
    --brief -y --wait

dx download ${INPUTS}/${COVAR}.csv

Rscript - <<-RSCRIPT
    library(tidyverse)
    extracted <- read_csv("${COVAR}.csv")
    covar_new <- extracted %>%
        rename_with(~str_replace(., "22009-0.", "PC")) %>%
        rename(array=\`22000-0.0\`, sex=\`31-0.0\`, age=\`21022-0.0\`) %>%
        mutate(array=ifelse(array < 0, 1, 2), FID=eid) %>%
        drop_na %>%
        select(FID, IID=eid, age, array, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
    write_tsv(covar_new, "${covar}")
RSCRIPT
