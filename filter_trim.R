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
    xlim(0, 2000)
ggsave("qualities.png", width=10, height=5)

cat("Filtering reads...\n")
if (!file.exists("filtered")) {
    filtered <- filterAndTrim(fastqs, file.path("filtered", basename(fastqs)),
                              trimLeft=10, truncLen=1500, maxN=2, truncQ=1,
                              rm.phix=TRUE, compress=TRUE, maxEE=100,
                              multithread=TRUE)
}
