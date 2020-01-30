library(readr)
source("R/plotting.R")
source("R/cleanup.R")
source("R/summarize.R")
source("R/reverse_operons.R")
#source("shiny/shinyfunctions.R")

####Data import####

all <- read_tsv("data/rawdata_tsv/old_txt/all_semiclean_20191219.txt")

lineages_map <- read_delim("data/acc_files/organisms2lineages_map_bae_20170828.txt",
                           delim="\t", col_names=T, trim_ws=T)

domains_rename <- read_delim("data/acc_files/domains_rename.txt",
                             delim="\t", col_names=TRUE)

domains_keep <- read_delim("data/acc_files/domains_keep.txt",
                           delim="\t", col_names=T)

query_domains <- read_delim("data/acc_files/query_domains.txt",
                            delim="\t", col_names=T)

### Cleanup will be removed once a cleaned file can be read ###

#
#
# # Cleanup Clusters
# all <- all %>%
#   cleanup_clust(repeat2s = TRUE, domains_keep, domains_rename)
#
# # Cleanup Species
# all <- all %>%
#   cleanup_species(remove_empty = FALSE)
#
# # Cleanup GenContext
# # Calls reverse_operons
# all <- all %>%
#   cleanup_gencontext(repeat2s = TRUE, remove_empty = FALSE, domains_rename)
#
# all <- all %>%
#   cleanup_domarch(repeat2s = TRUE, domains_rename)

colnames(all)[which(colnames(all) == "DomArch")] = "DomArch.ind"
colnames(all)[which(colnames(all) == "DomArch.orig")] = "DomArch.ind.orig"

colnames(all)[which(colnames(all) == "ClustName")] = "DomArch"

all <- all %>%
  cleanup_clust(domains_keep = domains_keep, domains_rename = domains_rename,repeat2s = FALSE)
colnames(all)[which(colnames(all) == "ClustName")] = "DomArch.repeats"

colnames(all)[which(colnames(all) == "ClustName.orig")] = "DomArch.orig"

all <- all %>% select(AccNum, DomArch, DomArch.repeats, DomArch.ind, GenContext, Lineage, Species, GeneName, Length, GCA_ID)

DUF1700 <- all %>% filter(grepl("DUF1700-ahelical",ignore.case = T, DomArch))
DUF1707 <- all %>% filter(grepl("DUF1707-SHOCT",ignore.case = T, DomArch))
pspa <- all%>% filter(grepl("pspa|snf7",ignore.case = T,DomArch))
pspb <- all%>% filter(grepl("pspb",ignore.case = T, DomArch))
pspc <- all%>% filter(grepl("pspc",ignore.case = T, DomArch))
pspm <- all%>% filter(grepl("pspm",ignore.case = T, DomArch))
pspn <- all%>% filter(grepl("pspn",ignore.case = T, DomArch))
liai_liaf = all%>% filter(grepl("LiaI-LiaF-TM",ignore.case = T, DomArch))
toast_rack = all%>%filter(grepl("Toast-rack",ignore.case = T, DomArch))
