# Align to 16S reference

library(mbtools)
library(ShortRead)
library(stringr)

files <- list.files("filtered", ".fastq.gz", full.names=T)

if (!file.exists("alignments")) {
    align_nanopore(files, "silva_species.fa.gz", threads=16)
}

alns <- list.files("alignments", ".bam", full.names=T)
counts <- count_nanopore(alns)
ref <- readFasta("silva_species.fa.gz")
ids <- as.character(id(ref))
annotations <- as.data.table(str_split_fixed(ids, " ", 4))
names(annotations) <- c("seqnames", "genus", "species", "strain")
counts_ann <- counts[annotations, on="seqnames"]

fwrite(counts_ann, counts.csv)
