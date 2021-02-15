# azimuth-references

Repo with workflows for generating azimuth reference objects. 

# Overview

This repository contains the Dockerfile and snakemake workflows that are used to generate the azimuth references that are hosted online. Each reference directory contains a `Snakefile` and associated scripts that can be run to regenerate each reference (and associated demo data) from publicly available download links of the underlying data. To run: 

```
snakemake --use-singularity --cores 1 all
```
