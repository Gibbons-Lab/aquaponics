# Nanopore full-length 16S sequencing benchmark

Nanopore long read sequencing is known to have large error rates. To this
end we compared Nanopore sequencing and Illumina sequencing on the same samples
from the biofilter of an established aquaponics system.

## Install

You will not need to install any additional software to what is outlined
in the [top repository](../README.md). All commands below are assume to be
run from this directory. If you are still in the top directory you can change this
with

```bash
cd nanopore_benchmark
```

## Denoise the Illumina sequencing reads and obtain ASVs

This is done with the `illumina.R` script for which you can adjust the
number of cores as before.

```bash
Rscript -e "options(mc.cores = 10); source('illumina.R')"
```

## Align and classify the nanopore reads

This is performed exactly the same way as the other samples in the manuscript
using the `nanopore.R` script.

```bash
Rscript -e "options(mc.cores = 10); source('nanopore.R')"
```

## Perform comparisons

This is again separated into a set of analysis steps which can be reproduced
in Rstudio by opening the mentioned notebooks.

- show improved performance of EM counting [notebook](methods.rmd) | [output](https://gibbons-lab.github.io/aquaponics/methods.nb.html)
- compare Nanopore and Illumina data [notebook](comparison.rmd) | [output](https://gibbons-lab.github.io/aquaponics/comparison.nb.html)
