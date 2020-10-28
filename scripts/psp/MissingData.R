# Samuel Chen
# Investigate which rows get dropped from the raw psp data

library(tidyverse)
source("R/summarize.R")

raw <- read_tsv("data/rawdata_tsv/all_raw.txt", col_names  = c("AccNum", "ClustID", "ClustName.orig", "GenContext.orig","DomArch.Pfam", "DomArch.orig", "-", "Length", "GeneName", "Lineage", "Species.orig", "GCA_ID", "Annotation", "GI"))

clean <- read_tsv("data/rawdata_tsv/all_clean.txt")
clean <- read_tsv("data/rawdata_tsv/all_clean1013.txt")

raw_acc <- raw$AccNum %>% unique()

clean_acc <- clean$AccNum %>% unique()

only_in_raw <- raw_acc[-which(raw_acc %in% clean_acc)]

length(only_in_raw)


only_in_raw_df <- raw[which(raw$AccNum %in% only_in_raw),]

# Frequencies of all Clustnames that are only contained in the raw DF
only_in_raw_CNs <- only_in_raw_df %>% group_by(ClustName.orig)%>%
  summarize(count = n()) %>% arrange(-count)

# Frequencies of all the domains
raw_wc <-only_in_raw_df %>%
  elements2words(column = "ClustName.orig", conversion_type = "da2doms") %>%
  words2wc()


# How many rows do each of the domains appear in?
raw_wc$Rows <-  map(frequent_raw_wc$words, function(x){
  only_in_raw_df %>% filter(grepl(x,ClustName.orig)) %>% nrow()
}) %>% unlist() %>% sort( decreasing = T)

view(only_in_raw_CNs)
view(raw_wc)

write_tsv(only_in_raw, "data/rawdata_tsv/DAs_only_in_raw.txt", col_names = T)
write_tsv(raw_wc, "data/rawdata_tsv/Words_only_in_raw.txt")