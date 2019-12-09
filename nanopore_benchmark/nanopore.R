library(mbtools)

pattern <- "(\\w+)\\.fastq"
annotations <- c("id")

files <- find_read_files("data/illumina.R", pattern = pattern,
                         annotations = annotations)
quals <- files %>% quality_control()
ggsave("figures/nanopore_quals.png", plot = quals$quality_plot, dpi = 300)

config <- list(
    preprocess = config_preprocess(
        trimLeft = 10,
        truncLen = 1500,
        maxEE = 200,
        out_dir = "nanopore_filtered",
        truncQ = 0,
        maxN = 2
    ),
    align = config_align(
        reference = "../data/silva_132_dna_nr99.fa.gz",
        alignment_dir = "alignments"
    ),
    count = config_count()
)


if (!file.exists("nanopore.rds")) {
    artifact <- quals %>% preprocess(config$preprocess) %>%
                          align_long_reads(config$align) %>%
                          count_references(config$count)
    saveRDS(artifact, "nanopore.rds")
} else {
    artifact <- readRDS("nanopore.rds")
}

if (!file.exists("nanopore_naive.rds")) {
    artifact_naive <- quals %>% preprocess(config$preprocess) %>%
        align_long_reads(config$align) %>%
        count_references(method = "naive")
    saveRDS(artifact_naive, "nanopore_naive.rds")
} else {
    artifact_naive <- readRDS("nanopore_naive.rds")
}

counts <- artifact$counts
counts[, method := "em"]
counts_naive <- artifact_naive$counts
counts_naive[, method := "naive"]
counts <- rbind(counts, counts_naive)
anns <- annotate_silva(config$align$reference)
counts_ann <- anns[counts, on = c(id = "reference")]
fwrite(counts_ann, "nanopore_counts.csv")

# Build the fully annotated phyloseq object
mat <- dcast(counts[method == "em"], reference ~ sample,
	     value.var = "counts", fill = 0)
ids <- mat$reference
mat[, reference := NULL]
mat <- as.matrix(mat)
rownames(mat) <- ids
setkey(counts_ann, "id")
taxa <- counts_ann[, .(id, kingdom, phylum, class, order,
                       family, genus, species)][, .SD[1], by = id]
setkey(taxa, "id")
taxa <- as.matrix(taxa[rownames(mat)])
rownames(taxa) <- taxa[, 1]
ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE),
               tax_table(taxa))
saveRDS(ps, "nanopore_phylo.rds")
