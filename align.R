# Align to 16S reference

library(mbtools)
library(ShortRead)
library(stringr)

files <- list.files("filtered", ".fastq.gz", full.names=T)

if (!file.exists("alignments")) {
    align_nanopore(files, "silva_dna_132.fa.gz", threads=16)
}

alns <- list.files("alignments", ".bam", full.names=T)
cat("Counting reads...\n")
counts <- count_transcripts(alns, "/proj/gibbons/refs/silva_dna_132.fna.gz")
accs <- fread("zcat data/taxmap_132.txt.gz",
              select=c("primaryAccession", "start", "stop", "taxid"))
accs[, "seqnames" := paste0(primaryAccession, ".", start, ".", stop)]
annotations <- fread("data/silva_taxonomy_132.csv")
annotations <- annotations[accs, on="taxid", allow.cartesian=TRUE]
counts_ann <- annotations[counts, on=c(seqnames="transcript")]
fwrite(counts_ann, "data/counts.csv")
rm(annotations)

# Build the fully annotated phyloseq object
mat <- dcast(counts, transcript ~ sample, value.var="counts", fill=0)
ids <- mat$transcript
mat[, transcript := NULL]
mat <- as.matrix(mat)
rownames(mat) <- ids
taxa <- dcast(unique(counts_ann[, .(seqnames, rank, name)]),
              seqnames ~ rank,
              value.var="name", fill=NA)
setkey(taxa, "seqnames")
taxa <- as.matrix(taxa[rownames(mat)])
rownames(taxa) <- taxa[, 1]
samples <- read.csv("barcodes.csv")
rownames(samples) <- samples$isb_id
ps <- phyloseq(otu_table(mat, taxa_are_rows=TRUE),
               tax_table(taxa),
               sample_data(samples))
saveRDS(ps, "data/taxonomy.rds")
