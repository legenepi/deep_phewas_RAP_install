#!/bin/bash

BASE=`dirname $0`

. ${BASE}/RAP.config

DOCKER_TAG=legenepi/deep_phewas
DOCKER_SAVE=DeepPheWAS.docker.tar.gz
EXTRA_OPTIONS=${BASE}/extraOptions.json
PLINK2_VERSION=avx2_20220514
PLINK2=plink2_linux_$PLINK2_VERSION

dx select $PROJECT_ID &&
dx mkdir -p $PROJECT_DIR &&
dx cd $PROJECT_DIR &&
dx ls | grep -w docker_build.sh && dx rm docker_build.sh
dx ls | grep -w Dockerfile && dx rm Dockerfile
dx ls | grep -w $DOCKER_SAVE && dx rm $DOCKER_SAVE
dx ls | grep -w fields-minimum.txt.gz && dx rm fields-minimum.txt.gz

dx upload ${BASE}/docker_build.sh
dx upload ${BASE}/Dockerfile 
dx ls | grep -w $PLINK2 || dx upload ${BASE}/$PLINK2 

dx run --brief -y --wait --watch swiss-army-knife \
	-iin=docker_build.sh \
	-iin=Dockerfile \
	-iin=plink2_linux_$PLINK2_VERSION \
	-icmd="cp docker_build.sh /tmp; chmod u+x /tmp/docker_build.sh; /tmp/docker_build.sh $DOCKER_TAG $PLINK2_VERSION $PACKAGE $DOCKER_SAVE" &&

DOCKER_FILE_ID=`dx ls --brief $DOCKER_SAVE` &&

cat <<-JSON > $EXTRA_OPTIONS
	{
	    "defaultRuntimeAttributes" : {
	        "docker" : "dx://${DOCKER_FILE_ID}"
	    }
	}	
JSON
