#####Generate Missing files
#####Samuel Chen
library(tidyverse)
#setwd("C:/Users/samue/Google_Drive/GitHub/the-approach")
source("R/reverse_operons.R")
source("R/summarize.R")
source("R/cleanup.R")

all_with_tax_gca <- read_tsv("data/rawdata_tsv/all_merged_gca_taxid.txt")

spec_orig <- select(all_with_tax_gca , AccNum ,Species.orig, Species.q) %>% mutate(species = Species.orig)
spec_new <- select(all_with_tax_gca, AccNum, Species.q,Species.orig) %>% mutate(species = Species.q)
colnames(spec_orig)[colnames(spec_orig)=="Species.orig"] <- "species"
colnames(spec_new)[colnames(spec_new) == "Species.q"   ] <- "species"

anti1 <- anti_join(spec_new, spec_orig , by = "species") # This has more
anti2 <- anti_join(spec_orig, spec_new , by = "species") # This has less

ant_spec <- anti1 %>% select(Species.orig,Species.q) %>% distinct()
ant_spec2 <- anti2 %>% select(Species.orig,Species.q) %>% distinct()


#write_tsv(ant_spec,"data/acc_files/prob_files/species_diff.txt")


gc_blank <- all_with_tax_gca %>% filter( grepl("^-$", GenContext.orig)) %>% distinct()
da_blank <- all_with_tax_gca %>% filter( grepl("^-$", DomArch.orig)) %>% distinct()

gca_blank <- all_with_tax_gca %>% filter( grepl("^-$", GCA_ID) | is.na(GCA_ID)) %>% distinct()
gca_non_euk <- gca_blank %>% filter(!grepl("^eukaryota", Lineage) )

lineage_fix <- all_with_tax_gca %>% group_by(Lineage) %>% summarise(count = n()) %>% arrange(+count)

#write_tsv(lineage_fix,"data/acc_files/prob_files/lineage_summary.txt")
#write_tsv(gc_blank, "data/acc_files/prob_files/gc_blank.txt")
#write_tsv(gca_blank, "data/acc_files/prob_files/gca_blank.txt")
#write_tsv(da_blank, "data/acc_files/prob_files/da_blank.txt")
