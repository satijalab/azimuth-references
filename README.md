# azimuth-references

Repo with workflows for generating azimuth reference objects. 

# Overview

This repository contains the Dockerfile and snakemake workflows that are used to generate the azimuth references that are hosted online. Each reference directory contains a `Snakefile` and associated scripts that can be run to regenerate each reference (and associated demo data) from publicly available download links of the underlying data. To run: 

```
snakemake --use-singularity --cores 1 all
```

# Reference format

The Azimuth package provides the `AzimuthReference` function to facilitate converting existing Seurat objects into the specific format expected by Azimuth. Details on the required reference format can be viewed [here](https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format). For examples starting with a Seurat object, see the `export.R` scripts in the workflows in this repo (e.g. [human pancreas](https://github.com/satijalab/azimuth-references/blob/master/human_pancreas/scripts/export.R)). 
