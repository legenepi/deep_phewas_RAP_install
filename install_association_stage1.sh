#!/bin/bash

java -jar dxCompiler-2.13.0.jar compile WDL/association_stage1.wdl -folder /deep_phewas/regenie_workflow -f -streamFiles perfile
