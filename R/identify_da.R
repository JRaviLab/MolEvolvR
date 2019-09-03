## Identifying domains and clusternames
library(tidyverse)
source("cleanup.R")

# Read master file
prot <- read_delim("data/rawdata_tsv/all_raw.txt",
                   delim="\t", col_names=T, comment="#", trim_ws=T,
                   escape_double=FALSE, na="NA")

# Identify unique cluster names
clustnames <- prot %>%
  select(ClustName) %>%
  distinct() %>% arrange() %>%
  repeat2s(by_column="ClustName")

write_tsv(clustnames, "data/acc_files/clustnames_uniq.txt")

# Identify unique domains to retain
domains_keep <- prot %>%
  select(Query) %>%
  distinct() %>% arrange()

write_tsv(domains_keep, "data/acc_files/domains_keep.txt")