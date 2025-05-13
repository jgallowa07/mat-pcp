#!/usr/bin/env bash

# This script is used to run the Python script with the specified arguments.
# these arguments were gathered from the output of the `dasm_tree` branch of s2trajectory pipeline.
pixi run python scripts/pcp_from_tree.py \
    data/results_GISAID-chron-5k-wo-mp-ref-2024-10-01/mat/filtered_tree.pb \
    data/results_GISAID-chron-5k-wo-mp-ref-2024-10-01/ref/ref.fa \
    data/results_GISAID-chron-5k-wo-mp-ref-2024-10-01/ref/edited_ref.gtf \
    results/all_pcps.csv
    
# data/results_GISAID-chron-5k-wo-mp-ref-2024-10-01/ref/coding_sites.csv \