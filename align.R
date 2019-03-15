# Align to 16S reference

library(mbtools)
library(stringr)

pattern <- "(\\w+).fastq.gz"
annotation <- c("id")

config <- list(
    preprocess = config_preprocess(
        trimLeft = 0,
        truncLen = 1500,
        maxEE = 200,
        out_dir = "data/filtered"
    ),
    align = config_align_long(
        reference = "silva_dna_132.fa.gz",
        threads = 20,
        alignment_dir = "data/alignments"
    ),
    count = config_count(
        threads = 20
    )
)

files <- find_read_files("data/raw", pattern = pattern, annotations = annotation)
quals <- quality_control(files)
ggsave("figures/qualities.png", plot = quals$quality_plot)

if (!file.exists("data/workflow.rds")) {
    artifact <- quals %>% preprocess(config$preprocess) %>%
                          align_long_reads(config$align) %>%
                          count_transcripts(config$count)
    saveRDS(artifact, "data/workflow.rds")
} else {
    artifact <- readRDS("data/workflow.rds")
}

counts <- artifact$counts
counts[, counts := round(counts)]
counts <- counts[counts > 0]
accs <- fread("zcat data/taxmap_132.txt.gz",
              select=c("primaryAccession", "start", "stop", "taxid"))
accs[, "seqnames" := paste0(primaryAccession, ".", start, ".", stop)]
annotations <- fread("data/silva_taxonomy_132.csv")
annotations <- annotations[accs, on = "taxid", allow.cartesian = TRUE]
counts_ann <- annotations[counts, on = c(seqnames = "transcript")]
fwrite(counts_ann, "data/counts.csv")
rm(annotations)

# Build the fully annotated phyloseq object
mat <- dcast(counts, transcript ~ sample, value.var = "counts", fill = 0)
ids <- mat$transcript
mat[, transcript := NULL]
mat <- as.matrix(mat)
rownames(mat) <- ids
taxa <- dcast(unique(counts_ann[, .(seqnames, rank, name)]),
              seqnames ~ rank,
              value.var = "name", fill = NA)
setkey(taxa, "seqnames")
taxa <- as.matrix(taxa[rownames(mat)])
rownames(taxa) <- taxa[, 1]
samples <- read.csv("data/barcodes.csv")
rownames(samples) <- samples$isb_id
ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE),
               tax_table(taxa),
               sample_data(samples))
saveRDS(ps, "data/taxonomy.rds")
