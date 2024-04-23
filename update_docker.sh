#!/bin/bash

BASE=`dirname $0`

. ${BASE}/RAP.config

DOCKER_TAG=legenepi/deep_phewas
DOCKER_SAVE=DeepPheWAS.docker.tar.gz
EXTRA_OPTIONS=${BASE}/extraOptions.json
PLINK2_VERSION=avx2_20220514
PLINK2=${BASE}/plink2_linux_$PLINK2_VERSION

dx select $PROJECT_ID &&
dx mkdir -p $PROJECT_DIR &&
dx cd $PROJECT_DIR &&
dx ls | grep -w Dockerfile && dx rm Dockerfile
dx ls | grep -w $DOCKER_SAVE && dx rm $DOCKER_SAVE

dx upload ${BASE}/Dockerfile &&
dx ls | grep -w $PLINK2 || dx upload $PLINK2 &&

dx run --brief -y --wait --watch swiss-army-knife \
	-iin=Dockerfile \
    -iin=plink2_linux_$PLINK2_VERSION \
	-icmd="docker build -t $DOCKER_TAG --build-arg PLINK2_VERSION=$PLINK2_VERSION --build-arg PACKAGE=$PACKAGE . && docker save $DOCKER_TAG | gzip > $DOCKER_SAVE" &&

DOCKER_FILE_ID=`dx ls --brief $DOCKER_SAVE` &&

cat <<-JSON > $EXTRA_OPTIONS
	{
	    "defaultRuntimeAttributes" : {
	        "docker" : "dx://${DOCKER_FILE_ID}"
	    }
	}	
JSON
