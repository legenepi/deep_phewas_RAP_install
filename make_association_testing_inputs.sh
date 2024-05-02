#!/bin/bash

BASE=`dirname $0`

. ${BASE}/RAP.config

${BASE}/make_association_testing_inputs.R $0
