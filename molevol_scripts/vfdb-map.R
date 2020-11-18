# Script to Identify Unique MAP VFs
# Created: Mar 22, 2019 | Janani Ravi  

##################
library(tidyverse)
library(here)
##################

myco_vfs_long <- read_tsv("data/mycobacteria_vfs_long.txt")
myco_vfs_wide <- read_tsv("data/mycobacteria_vfs_wide.txt")

relatedgene_count <- myco_vfs_long %>%
  filter(!is.na(gene_name)) %>%
  group_by(Related_genes) %>%
  summarise(RelGeneNum=n()) %>% arrange(-RelGeneNum)

myco_vfs_long %>%
  filter(!is.na(gene_name)) %>%
  count(gene_name, name="GeneNum")
myco_vfs_long %>%
  filter(!is.na(gene_name)) %>%
  count(Related_genes, name="RelGeneNum") %>% arrange(-RelGeneNum)
vf_counts <- myco_vfs_long %>%
  filter(!is.na(gene_name)) %>%
  count(Virulence_factors, sort=TRUE, name="VFNum") %>% arrange(VFNum)

## MAP genes w/ their prevalence across Mycobacterial species
mavium_vfs <- myco_vfs_long %>%
  filter(grepl("avium", x=species)) %>%  # Paratb
  filter(!is.na(gene_name)) %>%                     # Has genes
  left_join(relatedgene_count, by="Related_genes") %>%
  left_join(vf_counts, by="Virulence_factors") %>%
  arrange(RelGeneNum, Related_genes)
View(mavium_vfs)

## Top 25 non-mce genes in M. avium
mavium_vfs_nonmce <- mavium_vfs %>%
  filter(!grepl(pattern="mce", Related_genes)) %>%
  arrange(RelGeneNum, Related_genes) %>%
  top_n(n=-30, wt=RelGeneNum)
View(mavium_vfs_nonmce)

## Top 25 non-mce genes in other species
mavium_vfs_nonmce %>%
  select(Related_genes) %>%
  left_join(myco_vfs_long) %>%
  filter(!is.na(gene_name)) %>% View()

GeneNameQuery <- "rmt3"
myco_vfs_long %>%
  filter(!is.na(gene_name)) %>%
  filter(grepl(pattern=GeneNameQuery, Related_genes)) %>%
  select(species) %>% distinct()

here()
write_tsv(mavium_vfs_nonmce, path="data/mavium_vfs_nonmce.txt",col_names=T)
write_tsv(mavium_vfs, path="data/mavium_vfs.txt",col_names=T)
