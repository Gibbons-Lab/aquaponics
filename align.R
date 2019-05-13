# Align to 16S reference

library(mbtools)
library(stringr)
library(futile.logger)

pattern <- "(\\w+).fastq.gz"
annotation <- c("id")

flog.appender(appender.tee("align.log"))

config <- list(
    preprocess = config_preprocess(
        trimLeft = 10,
        truncLen = 1500,
        maxEE = 200,
        out_dir = "data/filtered",
        threads = 20,
        truncQ = 0,
        maxN = 2
    ),
    align = config_align(
        reference = "silva_132_dna_nr99.fa.gz",
        threads = 20,
        alignment_dir = "data/alignments"
    ),
    count = config_count(
        threads = 20
    )
)

files <- find_read_files("data/raw", pattern = pattern,
                         annotations = annotation)
quals <- quality_control(files)
ggsave("figures/qualities.png", plot = quals$quality_plot + xlim(0, 1800))

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
annotations <- annotate_silva(config$align$reference)
counts_ann <- annotations[counts, on = c(id = "transcript")]
fwrite(counts_ann, "data/counts.csv")
rm(annotations)

# Build the fully annotated phyloseq object
mat <- dcast(counts, transcript ~ sample, value.var = "counts", fill = 0)
ids <- mat$transcript
mat[, transcript := NULL]
mat <- as.matrix(mat)
rownames(mat) <- ids
taxa <- counts_ann[, .(id, kingdom, phylum, class, order, 
		       family, genus, species)][, .SD[1], by = id]
setkey(taxa, "id")
taxa <- as.matrix(taxa[rownames(mat)])
rownames(taxa) <- taxa[, 1]
samples <- read.csv("data/barcodes.csv")
rownames(samples) <- samples$isb_id
ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE),
               tax_table(taxa),
               sample_data(samples))
saveRDS(ps, "data/taxonomy.rds")
