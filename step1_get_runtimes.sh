#!/bin/bash

# This script was created by copilot

# Find all step1:main jobs
dx find jobs --name "step1_bt" -n 200 --json \
| jq -r '.[] | select(.function == "main") | .id' \
| while read -r jobid; do

    # Full JSON for this job
    json=$(dx describe "$jobid" --json)

    # --------------------------
    # Extract job name + job id
    # --------------------------
    jobname=$(echo "$json" | jq -r '.name')

    # --------------------------
    # Extract phenoColList + count items
    # --------------------------
    pheno=$(echo "$json" | jq -r '.input.phenoColList')

    if [ -n "$pheno" ] && [ "$pheno" != "null" ]; then
        num_pheno=$(echo "$pheno" | awk -F',' '{print NF}')
    else
        num_pheno=0
    fi

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # SKIP jobs where number of pheno entries = 0
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if [ "$num_pheno" -eq 0 ]; then
        continue
    fi

    # --------------------------
    # Extract runtime
    # --------------------------
    runtime=$(echo "$json" | jq -r '.runTime // empty')

    # If runTime not present, compute from timestamps
    if [ -z "$runtime" ] || [[ "$runtime" == "null" ]]; then
        start_ms=$(echo "$json" | jq -r '.startedRunning // empty')
        stop_ms=$(echo "$json" | jq -r '.stoppedRunning // empty')

        if [ -n "$start_ms" ] && [ -n "$stop_ms" ]; then
            start_s=$(( start_ms / 1000 ))
            stop_s=$(( stop_ms / 1000 ))
            diff=$(( stop_s - start_s ))
            runtime=$(date -u -d @"$diff" +%H:%M:%S)
        else
            runtime="NA"
        fi
    fi

    # --------------------------
    # Output table row
    # --------------------------
    echo -e "${num_pheno}\t${runtime}\t${jobname}\t${jobid}"

done
