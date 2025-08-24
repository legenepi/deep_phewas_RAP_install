#!/bin/bash

. RAP.config 

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
    --brief -y
