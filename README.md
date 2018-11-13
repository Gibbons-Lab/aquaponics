## :fish: :seedling: :tada: ISB Aquaponics

Analysis code and scripts for the Aquaponics collab with Jessica Day et. al.

## Preprocessing

`preprocess.py` just rearranges the reads so that all reads for one sample
get their own file. `dada2.R` was the original code to run the [DADA2](https://benjjneb.github.io/dada2/index.html)
analysis on the Nanopore reads, but that one did not give any sensible results.

We continued with normal trimming and filtering with the following parameters:

- left trim: 10bp
- length truncated to 1.5kbp
- max 2 Ns per sequence
- truncate when encountering first base with quality <1
- max. of 100 expected errors per seq (Illumina model)
- phiX removal

The reads that passed the filter (around 30 - 50%) were passed off to the alignment step.

## Mapping

Reads were mapped to >300k 16S reference sequences from the SILVA database (version 132). Mapping was performed with minimap2
using the Oxford Nanopore presets. Counting and resolution of multi-mapping was performed with an EM algorithm similar to what [kalliso uses](https://www.nature.com/articles/nbt.3519) assuming a single start position for each transcript. 

see https://github.com/Gibbons-Lab/mbtools for implementation details.


