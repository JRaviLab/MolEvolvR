## Purpose: To cleanup the PATRIC genome list based off GenBank/RefSeq
## Created: Aug 05, 2020
## Last modified: Oct 09, 2020
## Author(s): Janani Ravi

## Libraries
library(tidyverse)
library(here)

## Import data
genome_path <- here("data/genomes/")
# PATRIC
bac_patric <- read_csv(paste0(genome_path, "PATRIC-firmicutes-compl_refrep_good_public_202008.csv"),
                              # "PATRIC-bac-compl_refrep_good_public_202008.csv"),
                       col_names=T)
# Column rename | remove spaces
colnames_old <- colnames(bac_patric)
colnames_new <- str_replace_all(string=colnames_old,
                                pattern=" ",
                                replacement="_")
colnames(bac_patric) <- colnames_new

bac_patric_dis <- bac_patric %>%
  filter(!is.na(Disease))

# BacDive
human_bacdive <- read_csv2(paste0(genome_path, "Firmicutes_human_path_bacdive.csv"),
                       col_names=T, skip=8)
human_bacdive_spp <- human_bacdive$strains.Species %>% sort() %>% unique()

animal_bacdive <- read_csv2(paste0(genome_path, "Firmicutes_animal_path_bacdive.csv"),
                           col_names=T, skip=8)
animal_bacdive_spp <- animal_bacdive$strains.Species %>% sort() %>% unique()

bacdive_spp_pat <- paste0(paste0(animal_bacdive_spp, collapse="|"),
                          paste0(human_bacdive_spp, collapse="|"),
                          collapse="|")

bac_patric %>% filter(str_detect(pattern=bacdive_spp_pat,
                                 string=Genome_Name)) %>% select(Genome_Name)


# GenBank
colnames_assembly <- c("Assembly_Accession", "bioproject", "biosample", "wgs_master",
                      "refseq_category", "taxid", "species_taxid",
                      "organism_name", "infraspecific_name",  "isolate",
                      "version_status", "assembly_level", "release_type",
                      "genome_rep", "seq_rel_date", "asm_name","submitter",
                      "gbrs_paired_asm", "paired_asm_comp",
                      "ftp_path", "excluded_from_refseq",
                      "relation_to_type_material")
gb <- read_tsv(paste0(genome_path,
                          "assembly_summary_genbank.20190312.txt"),
                   col_names=colnames_assembly, skip=2)
# RefSeq
rs <- read_tsv(paste0(genome_path,
                          "assembly_summary_refseq.20190312.txt"),
                   col_names=colnames_assembly, skip=2)

bac_gb <- bac_patric %>%
  left_join(gb, by="Assembly_Accession")

