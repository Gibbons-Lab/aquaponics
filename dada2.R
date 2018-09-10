library(dada2)
library(data.table)
library(ggplot2)

samples <- fread("barcodes.csv")
names(samples) <- c("site", "flow_cell", "timepoint", "sample_no",
                    "barcode_no", "group", "id")
setkey(samples, id)
fastqs <- list.files("raw", pattern=".fastq.gz", full.names=T)
snames <- sub(".fastq.gz", "", basename(fastqs))
samples[snames, fastq := fastqs]

cat("Plotting qualities...\n")
pl <- plotQualityProfile(fastqs[sample(1:length(fastqs), 3)], n=10000) +
    xlim(0, 1800)
ggsave("qualities.png", width=10, height=5)

cat("Filtering reads...\n")
if (!file.exists("filtered")) {
    filtered <- filterAndTrim(fastqs, file.path("filtered", basename(fastqs)),
                              trimLeft=10, truncLen=1500, maxN=2, truncQ=1,
                              rm.phix=TRUE, compress=TRUE, maxEE=100,
                              multithread=TRUE)
}

cat("Estimating error models...\n")
if (!file.exists("error_model.rds")) {
    errors <- lapply(unique(samples[, flow_cell]), function(cell) {
        files <- file.path("filtered",
                           basename(samples[flow_cell == cell, fastq]))
        derep <- derepFastq(files)
        return(learnErrors(derep, multithread=T, BAND_SIZE=32))
    })
    names(errors) <- samples[, unique(flow_cell)]
    saveRDS(errors, "error_model.rds")
} else {
    errors <- readRDS("error_model.rds")
}

plots <- sapply(names(errors), function(cell) {
    plotErrors(errors[[cell]], nominalQ=T)
    ggsave(paste0("errors_", cell, "*.png"), width=15, height=10)
})

cat("Obtaining ASVs...\n")
dadas <- lapply(names(errors), function(cell) {
    files <- file.path("filtered", basename(samples[flow_cell == cell, fastq]))
    derep <- derepFastq(files)
    #dada(derep, err=errors[[cell]], multithread=TRUE)
})

seqtab <- mergeSequenceTables(tables=lapply(dadas, makeSequenceTable))



