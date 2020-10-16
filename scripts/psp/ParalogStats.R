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


all_num_lineages = all$Lineage %>% unique() %>% length()
all_num_species = all$Species %>% unique() %>% length()
queries = c("PspA", "Snf7", "PspB","PspC", "PspM", "PspN","DUF3046",
            "LiaI-LiaF-TM", "Toast_rack","ahelical", "SHOCT-bihelical")
lin_fracs = c()
species_fracs = c()

lin_fracs_all = c()
species_fracs_all = c()

# Do by both query_num_lineage
# and all_num_lineage

for(q in queries)
{
  print(paste0("-----",q,"-----"))
  query_data <- all %>% filter(grepl(q, DomArch))
  para <- find_paralogs(query_data)

  num_para_lins <- para$Lineage %>% unique() %>% length()

  query_num_lineages = query_data$Lineage %>% unique() %>% length()

  lin_frac <- num_para_lins/query_num_lineages

  lin_fracs <- append(lin_fracs, round(lin_frac,4) )

  lin_frac_all <- num_para_lins/all_num_lineages
  lin_fracs_all <- append(lin_fracs_all, round(lin_frac_all,4))

  print(paste0("Fraction of Lineages: ", round(lin_frac, 3)))

  # 2. In terms of species

  num_para_species <- para$Species %>% unique() %>% length()

  query_num_species = query_data$Species %>% unique() %>% length()

  species_frac <- num_para_species/query_num_species
  species_fracs <- append(species_fracs, round(species_frac,4) )

  species_frac_all <- num_para_species/all_num_species
  species_fracs_all <- append(species_fracs_all, round(species_frac_all,4) )


  print(paste0("Fraction of Species: ", round(species_frac, 3)))

}

para_stats <- data.frame("Query" = queries, "FracLinByQuery" = lin_fracs,
                         "FracSpeciesByQuery" = species_fracs,
                         "FracLinByAll" = lin_fracs_all,
                         "FracSpeciesByAll" = species_fracs_all
                         )



