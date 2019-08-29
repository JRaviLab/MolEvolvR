## Identifying domains and clusternames
library(tidyverse)
source("cleanup.R")


prot <- read_delim("data/rawdata_tsv/all_raw.txt",
                   delim="\t", col_names=T, comment="#", trim_ws=T,
                   escape_double=FALSE, na="NA")


# Identify unique cluster names
prot %>%
  select(ClustName) %>%
  distinct() %>% arrange() %>%
  repeat2s(by_column="ClustName") %>% # From ("cleanup.R")
  write_tsv("data/acc_files/clustnames_uniq.txt")

# Identify unique domains to retain
prot %>%
  select(Query) %>%
  distinct() %>% arrange() %>%
  write_tsv("data/acc_files/domains_keep.txt")