#####Generate Missing files
#####Samuel Chen; + JR has joined the party! :(
## issue: input file hasn't been merged with missing_* files before checking!!
## everything would have to be rerun soon!

library(tidyverse)
source("R/reverse_operons.R")
source("R/summarize.R")
source("R/cleanup.R")

all_with_tax_gca <- read_tsv("data/rawdata_tsv/all_merged_gca_taxid.txt")

###############################
## Species.orig vs Species.q ##
## Sam's take
spec_orig <- select(all_with_tax_gca , AccNum ,Species.orig, Species.q) %>%
  mutate(species = Species.orig)
spec_new <- select(all_with_tax_gca, AccNum, Species.q,Species.orig) %>%
  mutate(species = Species.q)
colnames(spec_orig)[colnames(spec_orig)=="Species.orig"] <- "species"
colnames(spec_new)[colnames(spec_new) == "Species.q"   ] <- "species"

anti1 <- anti_join(spec_new, spec_orig , by = "species") # This has more
anti2 <- anti_join(spec_orig, spec_new , by = "species") # This has less

ant_spec <- anti1 %>% select(Species.orig,Species.q) %>% distinct()
ant_spec2 <- anti2 %>% select(Species.orig,Species.q) %>% distinct()

#write_tsv(ant_spec,"data/acc_files/prob_files/species_diff.txt")

## J's take on Species.orig vs Species.q
species.q_diff <- prot %>% filter(!Species.q %in% Species.orig)
species.orig_diff <- prot %>% filter(!Species.orig %in% Species.q)

Spp.q.na <- prot %>% filter(is.na(Species.q))
Spp.orig.na <- prot %>% filter(is.na(Species.orig))

# NAs in Species.q
# 163
sum(is.na(prot$Species.q))
# Species.orig values for which Species.q has nothing
prot$Species.orig[is.na(prot$Species.q)] %>% sort() %>% unique()

prot.sp.q.na <- prot %>%
  filter(is.na(Species.q))

# Last written on: Oct 10, 2019
#write_tsv(x=prot.sp.q.na, path="data/acc_files/prob_files/species.q.na.txt")

# Species.q with no caps
Sp.q.small <- prot %>%
  filter(!grepl(pattern="[A-Z]", x=Species.q)) %>%
  group_by(Species.q) %>%
  arrange() %>%
  summarise(tot=n())
# Sp.q.small: this is OK except for the first one -- bacterium ...

# Sp.q.small <- prot %>%
#   filter(!grepl(pattern="[A-Z]", x=Species.q)) %>%
#   write_tsv("data/acc_files/prob_files/species.q.small_na.txt")

# Species.q with CAPS -- checking...
Sp.q.caps <- prot %>%
  filter(grepl(pattern="[A-Z]", x=Species.q)) %>%
  group_by(Species.q) %>%
  arrange() %>%
  summarise(tot=n())
# Terrible: there are ones like "African clawed frog" -- so no luck!

#########################################
## Empty rows in GC, DA, GCA, GCA_non_euk
## NEED TO INCLUDE NA and "" as well
empty_match <- c("^$", "^ $", "^-$")
gc_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), GenContext.orig) |
           is.na(GenContext.orig))

da_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), DomArch.orig) |
           is.na(DomArch.orig))

gca_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), GCA_ID) |
           is.na(GCA_ID))
gca_non_euk <- gca_missing %>% filter(!grepl("^eukaryota", Lineage))

taxid_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), TaxID) |
           is.na(TaxID))

spp_oldnew_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), Species.orig) |
           is.na(Species.orig) |
           grepl(paste(empty_match,collapse="|"), Species.q) |
           is.na(Species.q))

lin_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), Lineage) |
           is.na(Lineage))

gi_missing <- all_with_tax_gca %>%
  filter(grepl(paste(empty_match,collapse="|"), GI) |
           is.na(GI))

# Summarizing by lineage
# ?? wonder why??
lineage_fix <- all_with_tax_gca %>% group_by(Lineage) %>% summarise(count = n()) %>% arrange(+count)

# Last written: Oct 11, 2019
#write_tsv(lineage_fix,"data/acc_files/prob_files/lineage_summary.txt")
# write_tsv(gc_missing, "data/acc_files/prob_files/gc_missing.txt")
# write_tsv(gca_missing, "data/acc_files/prob_files/gca_missing.txt")
# write_tsv(da_missing, "data/acc_files/prob_files/da_missing.txt")
# write_tsv(taxid_missing, "data/acc_files/prob_files/taxid_missing.txt")
# write_tsv(spp_oldnew_missing,
#           "data/acc_files/prob_files/spp_oldnew_missing.txt")
# write_tsv(gi_missing, "data/acc_files/prob_files/gi_missing.txt")
