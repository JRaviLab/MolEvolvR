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
GCA2Lins <- function(prot_data, assembly_path = "data/acc_files/assembly_summary20201018.txt",
                       lineagelookup_path = "data/lineagelookup.txt" )
{
  #' Function that maps GCA_ID to taxid, and that taxid to a lineage
  #' Note: currently configured to have at most kingdom and phylum
  #' @author Samuel Chen
  #' @param prot_data Dataframe containing a column `GCA_ID`
  #' @param assembly_path String of the path to the assembly_summary path
  #' This file can be generated using the "DownloadAssemblySummary()" function
  #' @param lineagelookup_path String of the path to the lineage lookup file
  #' (taxid to lineage mapping). This file can be generated using the
  #' "create_lineage_lookup()" function

  assembly_summary <- fread(assembly_path ,sep = "\t")
  assembly_summary <- setnames(assembly_summary, "assembly_accession", "GCA_ID")

  mergedTax <- merge(x = prot_data,y = assembly_summary,by = "GCA_ID", all.x = T)

  lineage_map <- fread(lineagelookup_path, sep = "\t")

  mergedLins <- merge(mergedTax, lineage_map, by.x = "taxid", by.y="tax_id",
                      all.x = T)

  return(mergedLins)
}



