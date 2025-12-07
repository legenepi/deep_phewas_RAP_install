#!/bin/bash

. RAP.config

if [ ! -s $EXTRA_OPTIONS ]; then
    echo "$EXTRA_OPTIONS not found, running update_docker.sh"
    ./update_docker.sh
fi

[ -s $DXCOMPILER ] || (echo "Downloading dxCompiler"; wget $DXCOMPILER_URL -O $DXCOMPILER)

echo "Compiling WDL"
for i in $WDL; do
    java -jar $DXCOMPILER compile WDL/$i -extras $EXTRA_OPTIONS -project $PROJECT_ID -folder $PROJECT_DIR -streamFiles perfile -f
done
