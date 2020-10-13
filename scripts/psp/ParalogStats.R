# Samuel Chen
# Find the fraction of genomes that contain paralogs based in terms of:
# 1) lineages
# 2) species

library(tidyverse)
source("R/summarize.R")
# Paralogs: instances where the there are more than 1 gene in the genome: ie)

# Fraction of genomes that contain paralogs in terms of lineages and species
# 1. In terms of lineages

all <- read_tsv("data/rawdata_tsv/all_clean1013.txt")

para <- find_paralogs(all)

num_para_lins <- para$Lineage %>% unique() %>% length()

num_lineages = all$Lineage %>% unique() %>% length()

lin_frac <- num_para_lins/num_lineages

print(paste0("Fraction of Lineages: ", round(lin_frac, 3)))

# 2. In terms of species

num_para_species <- para$Species %>% unique() %>% length()

num_species = all$Species %>% unique() %>% length()

species_frac <- num_para_species/num_species

print(paste0("Fraction of Species: ", round(species_frac, 3)))

