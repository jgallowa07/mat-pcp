# `mat-pcp`

A simple pipeline for filtering, optimizing, and extraction pcp's from a MAT tree.

## Running the pipeline 

You must have an environment with `snakemake`, and `conda` available. Additionally, the config is setup to read from a
MAT data we don't store here. See [config.yaml](config/config.yaml) for more.
```
snakemake -j 4 --use-conda
```