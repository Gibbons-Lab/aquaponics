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
counts_ann <- annotations[counts, on="seqnames"]

fwrite(counts_ann, "counts.csv")

# Build the fully annotated phyloseq object
mat <- dcast(counts, seqnames ~ sample, value.var="counts", fill=0)
seqnames <- mat$seqnames
mat[, seqnames := NULL]
mat <- as.matrix(mat)
rownames(mat) <- seqnames
setkey(annotations, seqnames)
taxa <- as.matrix(annotations[rownames(mat)])
rownames(taxa) <- taxa[, 1]
samples <- read.csv("barcodes.csv")
rownames(samples) <- samples$isb_id
ps <- phyloseq(otu_table(mat, taxa_are_rows=TRUE),
               tax_table(taxa)
               sample_data(samples))
saveRDS(ps, "taxonomy.rds")
