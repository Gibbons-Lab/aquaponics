library(mbtools)
library(Biostrings)

pattern <- "(\\w+_\\d+)_(\\w+)_R(\\d+)\\.fastq"
annotations <- c("id", "barcode", "direction")

files <- find_read_files("data/illumina", pattern = pattern,
                         annotations = annotations)
quals <- files %>% quality_control()
ggsave("illumina_quals.png", plot = quals$quality_plot, dpi = 300)

conf <- list(
    filter = config_preprocess(
        trimLeft = 10,
        truncLen = c(230, 230),
        out_dir = "data/illumina_filtered"
    ),
    denoise = config_denoise()
)

if (!file.exists("illumina.rds")) {
    denoised <- quals %>% preprocess(conf$filter) %>%
       	        denoise(conf$denoise) %>% saveRDS("illumina.rds")
}
denoised <- readRDS("illumina.rds")

saveRDS(as_phyloseq(denoised), "illumina_phylo.rds")
