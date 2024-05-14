suppressPackageStartupMessages({
library(readr, quietly = T)
library(dplyr, quietly = T)
})

args = commandArgs(trailingOnly = TRUE)

workdir <- args[1]
threshold <- as.numeric(args[2])
sample_id <- args[3]
chrs <- args[4]
starts <- as.numeric(args[5])
ends <- as.numeric(args[6])
prefix <- args[7]

cov <- read.table(gzfile(file.path(workdir, paste0(sample_id, ".sam.coverage.gz"))), stringsAsFactors = F)
colnames(cov) <- c("chr", "pos", "reads")
cov_filt <- cov %>%
    filter(chr == chrs) %>%
    filter(pos >= starts & pos <= ends) %>%
    mutate(reads_count = n()) %>%
    filter(reads >= threshold) %>%
    mutate(reads_pass = n(),
           reads_perc = round(reads_pass / reads_count * 100, 2),
           reads_mean = mean(reads, na.rm = T),
           reads_median = median(reads, na.rm = T)) %>%
    select(-c(chr, pos, reads)) %>% distinct() %>%
mutate(chr = paste0(chrs, ":", starts, "-", ends))

cov_filt$sample_id <- sample_id

write.table(cov_filt, paste0(sample_id, "_", prefix, "_th_", threshold, "X_cov_summary.txt"), col.names=F, row.names=F, quote=F, sep="\t")