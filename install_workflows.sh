#!/bin/bash

BASE=`dirname $0`

. ${BASE}/RAP.config

EXTRA_OPTIONS=${BASE}/extraOptions.json
DXCOMPILER=/tmp/dxCompiler.jar
WDL="phenotype_generation.wdl association_testing.wdl"

[ -s $DXCOMPILER ] || wget https://github.com/dnanexus/dxCompiler/releases/download/2.11.4/dxCompiler-2.11.4.jar -O $DXCOMPILER &&
for i in $WDL; do
    java -jar $DXCOMPILER compile ${BASE}/WDL/$i -extras $EXTRA_OPTIONS -project $PROJECT_ID -folder $PROJECT_DIR -streamFiles perfile -f
done
