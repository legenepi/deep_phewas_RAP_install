#!/bin/bash

. RAP.config

dx cd $PROJECT_DIR 

if ! dx ls $FIELDS_MINIMUM; then
    echo "You need to run update_docker.sh first to extract $FIELDS_MINIMUM"
    exit 1
fi

dx mkdir -p $INPUTS &&
    (dx ls ${INPUTS}/$FIELDS_USE && dx rm -a ${INPUTS}/$FIELDS_USE) &&
    dx extract_dataset "$DATASET" --list-fields | cut -f1 | cut -f2 -d. |
      grep -Ef <(dx cat $FIELDS_MINIMUM | zcat | awk 'NR == 1 { print "eid"} {print "p"$0"($|_)"}') > $FIELDS_USE &&
    dx upload --path ${INPUTS}/ $FIELDS_USE &&

dx run table-exporter \
    -idataset_or_cohort_or_dashboard=$DATASET \
    -ientity="participant" \
    -ifield_names_file_txt="${INPUTS}/$FIELDS_USE" \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ioutput="$MINIMUM_DATA" \
    --destination $INPUTS \
    --instance-type mem2_ssd2_x16 \
    --name extract_fields \
    --brief -y

for ENTITY in $ENTITIES; do
    dx run table-exporter \
        -idataset_or_cohort_or_dashboard=$DATASET \
        -ientity="$ENTITY" \
        -ioutput_format=TSV \
        -icoding_option=RAW \
        -iheader_style=UKB-FORMAT \
        -ioutput="$ENTITY" \
        --destination $INPUTS \
        --name extract_$ENTITY \
        --brief -y
done
