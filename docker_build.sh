#!/bin/bash

DOCKER_TAG=$1
PLINK2_VERSION=$2
PACKAGE=$3
DOCKER_SAVE=$4

docker build -t $DOCKER_TAG --build-arg PLINK2_VERSION=$PLINK2_VERSION --build-arg PACKAGE=$PACKAGE . 
docker save $DOCKER_TAG | gzip > $DOCKER_SAVE 
FILE=`docker run $DOCKER_TAG Rscript '-e' 'cat(system.file("extdata", "fields-minimum.txt.gz", package = "DeepPheWAS"))'`
CONTAINER_ID=`docker run -d $DOCKER_TAG`
docker cp $CONTAINER_ID:$FILE .
docker stop $CONTAINER_ID
