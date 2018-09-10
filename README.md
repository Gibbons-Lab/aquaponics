## :fish: :seedling: :tada: ISB Aquaponics

Analysis code and scripts for the Aquaponics collab with Jessica Day et. al.

## Preprocessing

`preprocess.py` just rearranges the reads so that all reads for one sample
get their own file. `dada2.R` was the original code to run the [DADA2](https://benjjneb.github.io/dada2/index.html)
analysis on the Nanopore reads, but that one did not give any sensible results.
