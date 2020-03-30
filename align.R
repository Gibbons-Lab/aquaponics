# Align to 16S reference

library(mbtools)
library(stringr)
library(futile.logger)

pattern <- "(\\w+).fastq.gz"
annotation <- c("id")

file.remove("align.log")
flog.appender(appender.tee("align.log"))

config <- list(
    preprocess = config_preprocess(
        trimLeft = 10,
        truncLen = 1500,
        maxEE = 200,
        out_dir = "data/filtered",
        truncQ = 0,
        maxN = 2
    ),
    chimera = config_chimera(
        out_dir = "data/non_chimeric"
    ),
    align = config_align(
        reference = "data/silva_132_dna_nr99.fa.gz",
        alignment_dir = "data/alignments",
        limited_memory = FALSE
    ),
    count = config_count(
        weights = TRUE
    )
)

if (!file.exists("data/workflow.rds")) {
    files <- find_read_files("data/raw", pattern = pattern,
                             annotations = annotation)
    quals <- quality_control(files)
    ggsave("figures/qualities.png", plot = quals$quality_plot + xlim(0, 1800))
    artifact <- quals %>% preprocess(config$preprocess) %>%
                          remove_chimeras(config$chimera) %>%
                          align_long_reads(config$align) %>%
                          count_references(config$count)
    saveRDS(artifact, "data/workflow.rds")
} else {
    artifact <- readRDS("data/workflow.rds")
}

counts <- artifact$counts
counts[, counts := round(counts)]
counts <- counts[counts > 0 & !is.na(reference)]
annotations <- annotate_silva(config$align$reference)
counts_ann <- annotations[counts, on = c(id = "reference")]
fwrite(counts_ann, "data/counts.csv")
rm(annotations)

# Build the fully annotated phyloseq object
mat <- dcast(counts, reference ~ sample, value.var = "counts", fill = 0)
ids <- mat$reference
mat[, reference := NULL]
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
