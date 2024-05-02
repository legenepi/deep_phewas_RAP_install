#!/bin/bash

. RAP.config

TAB_DATA="death.tsv
death_cause.tsv
fields_use.txt
gp_clinical.tsv
gp_scripts.tsv
hesin.tsv
hesin_diag.tsv
hesin_oper.tsv
minimum_tab_data.csv"

for i in $TAB_DATA; do
    if ! dx ls ${INPUTS}/$i > /dev/null; then
        echo "Required input $i missing from ${INPUTS}, have you run extract_fields.sh?"
        exit 1
    fi
done

export PROJECT_ID INPUTS

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER

./R/`basename ${0/sh/R}` \
    --data_files "${INPUTS}/minimum_tab_data.csv" \
    --GPC "${INPUTS}/gp_clinical.tsv" \
    --GPP "${INPUTS}/gp_scripts.tsv" \
    --hesin_diag "${INPUTS}/hesin_diag.tsv" \
    --HESIN "${INPUTS}/hesin.tsv" \
    --hesin_oper "${INPUTS}/hesin_oper.tsv" \
    --death_cause "${INPUTS}/death_cause.tsv" \
    --death "${INPUTS}/death.tsv" \
    --king_coef "$KING_COEF" \
    --out "${PHENOTYPES}.json" \
    $* &&
    grep -q -- --help <<<"$*" ||
    java -jar $DXCOMPILER compile WDL/phenotype_generation.wdl -compileMode IR -inputs ${PHENOTYPES}.json
