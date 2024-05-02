#!/bin/bash

. RAP.config

WDL="phenotype_generation.wdl association_testing.wdl"

if [ ! -s $EXTRA_OPTIONS ]; then
    echo "$EXTRA_OPTIONS not found, have you run update_docker.sh?"
    exit 1
fi

[ -s $DXCOMPILER ] || wget $DXCOMPILER_URL -O $DXCOMPILER    

for i in $WDL; do
    java -jar $DXCOMPILER compile WDL/$i -extras $EXTRA_OPTIONS -project $PROJECT_ID -folder $PROJECT_DIR -streamFiles perfile -f
done
