#!/bin/bash

. RAP.config

dx cd $PROJECT_DIR 

if ! dx ls $FIELDS_MINIMUM; then
    echo "You need to run update_docker.sh first to extract $FIELDS_MINIMUM"
    exit 1
fi


if ! dx ls $DATASET; then
    echo "$DATASET not found, have you specified the correct dataset in options.config?"
    exit 1
fi

dx mkdir -p $INPUTS &&
    dx ls ${INPUTS}/$FIELDS_USE 2> /dev/null && dx rm -a ${INPUTS}/$FIELDS_USE

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
    --ignore-reuse \
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
        --ignore-reuse \
        --brief -y
done
