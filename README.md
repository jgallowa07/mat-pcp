# `mat-pcp`

A simple pipeline for filtering, optimizing, and extraction pcp's from a MAT tree. Notably, the pipeline allows you to select a datetime cutoff such that only samples before a given date are included. After filtering of samples the pipeline, it re-optimizes such that we have a tree of maximum parsimony. 

**NOTE** Currently, this pipeline is confined to the SCV2 tree and extracts the spike gene - but is easy to generalize if/when we're ready.

## Running the pipeline 

You must have an environment with `snakemake`, and `conda` available. Additionally, the config is setup to read from a
MAT data we don't store here. See [config.yaml](config/config.yaml) for more.
```
snakemake -j 4 --use-conda
```