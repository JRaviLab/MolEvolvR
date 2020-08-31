library(tidyverse)
library(biomartr)

DownloadAssemblySummary <- function(outpath, keep = c("assembly_accession", "taxid", "species_taxid", "organism_name"))
{
  #' Download the combined assembly summaries of genbank and refseq
  #' @author Samuel Chen
  #' @param outpath String of path where the assembly summary file should be written
  #' @param keep Character vector containing which columns should be retained and downloaded

  assembly_kingdom_genbank <- getKingdomAssemblySummary("genbank")
  assembly_kingdom_refseq <- getKingdomAssemblySummary("refseq")

  if(keep == "all")
  {
    assembly_all <- bind_rows(assembly_kingdom_genbank,assembly_kingdom_refseq)
  }
  else
  {
    assembly_all <- bind_rows(assembly_kingdom_genbank,assembly_kingdom_refseq) %>%
      select(all_of(keep))
  }

  write_tsv(assembly_all, outpath)
}




# Go from the GCA_ID column to tax IDs using the assembly file
GCA2Lins <- function(prot_data, assembly_path = "data/acc_files/assembly_summary20200706.txt", lineage_path = "data/acc_files/rankedlineage.dmp" )
{
  #' Function that maps GCA_ID to taxid, and that taxid to a lineage
  #' Note: currently configured to have at most kingdom and phylum
  #' @author Samuel Chen
  #' @param prot_data Dataframe containing a column `GCA_ID`
  #' @param assembly_path String of the path to the assembly_summary path
  #' @param lineage_path String of the path to the rankedLineage dump file (taxid to lineage mapping)

  shorten_NA <- function(Lineage)
  {
    first_NA = str_locate(Lineage, "NA")[1]
    if(is.na(first_NA) )
    {
      # No NAs
      # print(Lineage)
      return(Lineage)
    }
    else
    {
      shortened = substr(Lineage,1,(first_NA-1))
      return(shortened)
    }
  }

  assembly_summary <- read_tsv("data/acc_files/assembly_summary20200706.txt")

  assembly_summary <- rename(assembly_summary, "GCA_ID" = "assembly_accession")

  mergedTax <- prot_data %>% left_join(assembly_summary,by = "GCA_ID")

  dont_include = c("X2", "X4", "X6", "X8","X10","X12", "X14","X16","X17", "X18","X20")
  rankedLins <- read_tsv("data/acc_files/rankedlineage.dmp", col_names = F) %>%
    select(-all_of(dont_include)) %>% rename(taxid = X1)

  merged_lins <- left_join(mergedTax, rankedLins, by = "taxid")

  merged_lins <- merged_lins %>% unite(col = 'Lineage', X19,X15,#X19:X5
                                       sep = "_") %>%
    mutate(Lineage = map(Lineage,shorten_NA)) %>%
    unite(Lineage, Lineage, X3 )

  return(merged_lins)
}



