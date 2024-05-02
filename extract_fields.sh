#!/bin/bash

BASE=`dirname $0`

. ${BASE}/RAP.config

DEST=${PROJECT_DIR}/inputs

dx cd $PROJECT_DIR 

if ! dx ls fields-minimum.txt.gz; then
    echo "You need to run update_docker.sh first to extract fields-minimum.txt.gz"
    exit 1
fi

dx mkdir -p $DEST &&
    (dx ls ${DEST}/fields_use.txt && dx rm -a ${DEST}/fields_use.txt) &&
    dx extract_dataset "$DATASET" --list-fields | cut -f1 | cut -f2 -d. |
      grep -Ef <(dx cat fields-minimum.txt.gz | zcat | awk 'NR == 1 { print "eid"} {print "p"$0"($|_)"}') > fields_use.txt &&
    dx upload --path ${DEST}/ fields_use.txt &&

dx run table-exporter \
    -idataset_or_cohort_or_dashboard=$DATASET \
    -ientity="participant" \
    -ifield_names_file_txt="${DEST}/fields_use.txt" \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ioutput="minimum_tab_data" \
    --destination $DEST \
    --instance-type mem2_ssd2_x16 \
    --name extract_fields \
    --brief -y

for ENTITY in death death_cause gp_clinical gp_scripts hesin hesin_diag hesin_oper; do
    dx run table-exporter \
        -idataset_or_cohort_or_dashboard=$DATASET \
        -ientity="$ENTITY" \
        -ioutput_format=TSV \
        -icoding_option=RAW \
        -iheader_style=UKB-FORMAT \
        -ioutput="$ENTITY" \
        --destination $DEST \
        --name extract_$ENTITY \
        --brief -y
done
