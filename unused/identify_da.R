## Identifying domains and clusternames
library(tidyverse); library(wordcloud)
source("R/cleanup.R")

# Read master file
prot <- read_delim("data/rawdata_tsv/all_raw.txt",
                   delim="\t", col_names=T, comment="#", trim_ws=T,
                   escape_double=FALSE, na="NA")

# Read domains to keep
domains_keep <- read_delim("data/acc_files/domains_keep.txt",
                   delim="\t", col_names=T, comment="#", trim_ws=T,
                   escape_double=FALSE, na="NA")

# Identify unique domains to retain
# Not used after Aug 2019
# domains_keep <- prot %>%
#   select(Query) %>%
#   distinct() %>% arrange()

# write_tsv(domains_keep, "data/acc_files/domains_keep.txt")


# Identify unique cluster names
clustnames <- prot %>%
  select(ClustName) %>%
  arrange() %>% #distinct() %>%
  repeat2s(by_column="ClustName") %>%
  count(ClustName, sort=T)

clustnames.keep <- clustnames %>%
  filter(grepl(x=ClustName,
               pattern=paste(domains_keep$domains,
                             collapse="|")))
clustnames.rm <- clustnames %>%
  filter(!grepl(x=ClustName,
                pattern=paste(domains_keep$domains,
                              collapse="|")))

# Last written: Sep 03, 2019
write_tsv(clustnames.keep, "data/acc_files/clustnames_keep.txt")
write_tsv(clustnames.rm, "data/acc_files/clustnames_remove.txt")

## Wordcloud by query
## Temp frequency cutoff: 10
clust.plot <- clustnames.keep %>%
  # filter(grepl(pattern="PspA|Snf7", x=ClustName))
  # filter(grepl(pattern="Toast-rack", x=ClustName))
  # filter(grepl(pattern="PspC|PspB", x=ClustName))
  filter(grepl(pattern="DUF170", x=ClustName))
  # filter(grepl(pattern="DUF3046|PspN|PspM", x=ClustName))

wordcloud(words=clust.plot$ClustName, freq=clust.plot$n,
          min.freq=10, #rot.per=0.5,
          colors=(RColorBrewer::brewer.pal(n=9, "Dark2")),
          scale=c(3,0.5))
