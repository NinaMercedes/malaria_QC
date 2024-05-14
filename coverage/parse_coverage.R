suppressPackageStartupMessages({
library(readr, quietly = T)
library(dplyr, quietly = T)
library(data.table, quietly = T)
})


threshold <- 5

coverage <- list()

temp <- list.files (".", full.names = FALSE, pattern = "\\.samtools_depth.csv$")

for (file_name in 1:length(temp)){
  cov <- read.table(temp[file_name])
  sample_id <- gsub(".bqsr.cram.samtools_depth.csv","", file_name)
  colnames(cov) <- c("chr", "pos", "reads")
  cov_filt <- cov %>% filter(grepl("^P", chr)) %>%
      mutate(reads_count = n()) %>%
      filter(reads >= threshold) %>%
      mutate(reads_pass = n(),
             reads_perc = round(reads_pass/reads_count*100,2),
             reads_mean = round(mean(reads, na.rm = T),3),
             reads_median = median(reads, na.rm = T)) %>%
      select(-c(chr, pos, reads)) %>% distinct() %>%
  mutate(chr = "overall")
  
  cov_filt_chr <- cov %>% filter(grepl("^P", chr)) %>%
      group_by(chr) %>%
      mutate(reads_count = n()) %>%
      filter(reads >= threshold) %>%
      mutate(reads_pass = n(),
             reads_perc = round(reads_pass/reads_count*100,2),
             reads_mean = round(mean(reads, na.rm = T),3),
             reads_median = median(reads, na.rm = T)) %>%
      select(-c(pos, reads)) %>% distinct()
  
  combined_cov <-  cov_filt_chr %>% full_join(cov_filt)
  combined_cov$sample_id <- sample_id
  coverage[[file_name]] <- combined_cov
}
out <- rbindlist(coverage, fill = TRUE)

write.table(out, "coverage_all.txt", col.names=T, row.names=F, quote=F, sep="\t")
