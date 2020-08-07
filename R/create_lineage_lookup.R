library(tidyverse)
library(here)
# library(biomartr)


create_lineage_lookup <- function(lineage_file = here("data/rankedlineage.dmp"),
                                  outfile, taxonomic_rank = "phylum")
{
  #' Create a look up table that goes from TaxID, to Lineage
  #' @author Samuel Chen
  #' @param lineage_file Path to the rankedlineage.dmp file containing taxid's and their
  #' corresponding taxonomic rank. rankedlineage.dmp can be downloaded at
  #' https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
  #'@param outfile File the resulting lineage lookup table should be written to
  #'@param taxonomic_rank The upperbound of taxonomic rank that the lineage includes. The lineaege will
  #'include superkingdom>...>taxonomic_rank.

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
      # Only NA's
      if(first_NA == 1)
      {
        shortened = ""
      }
      else
      {
        shortened = substr(Lineage,1,(first_NA-2))
      }
      return(shortened)
    }
  }

  rankedlins_cols <- c("tax_id", "tax_name", "species", "genus",
                       "family", "order", "class", "phylum", "kingdom", "superkingdom", "NA")

  lineage_file <- here("data/rankedlineage.dmp")
  rankedLins <- read_file(lineage_file) %>%
    str_replace_all(pattern = "\\t\\|\\t|\\t\\|", "\t") %>%
    read_tsv( col_names = rankedlins_cols) %>%
    select(-"NA")

  rankedlins_cols_nona <- c("tax_id", "tax_name", "species", "genus",
                            "family", "order", "class", "phylum", "kingdom", "superkingdom")

  # kingdom column is empty?
  # rankedLins %>% group_by(kingdom) %>% summarize(n())

  taxonomy_ranks <- c("superkingdom",# "kingdom",
                         "phylum", "class", "order", "family", "genus", "species")

  tax_rank_index <- switch( taxonomic_rank,
    "superkingdom" = 1,
    # "kindom" = 2,
    "phylum" = 2,
    "class" = 3,
    "order" = 4,
    "family" = 5,
    "genus" = 6,
    "species" = 7
    )

  combined_taxonomy <- taxonomy_ranks[1:tax_rank_index]

  if(tax_rank_index < length(taxonomy_ranks))
  {
    rankedLins <- rankedLins %>%
      select(- all_of(taxonomy_ranks[ (tax_rank_index+1) :length(taxonomy_ranks)]), -kingdom)
  }
  else
  {
    rankedLins <- rankedLins %>%
      select(-kingdom)
  }


  # Takes a while (2million rows after all)
  rankedLinsCombined <- rankedLins %>%
    unite(col = 'Lineage', all_of(combined_taxonomy), sep = ">") %>%
    mutate(Lineage = map(Lineage, shorten_NA))

  write_tsv(rankedLinsCombined, outfile)
}



#' CreateLineageLookup <- function(assembly_path, updateAssembly = FALSE, file_type = "tsv")
#' {
#'   #' Create a look up table that goes from GCA_ID, to TaxID, to Lineage
#'   #' @author Samuel Chen
#'   #' @param assembly_path Path to the assembly summary to be read, and/or where
#'   #' the generated assembly file is stored based on the 'updateAssembly' parameter
#'   #' @param updateAssembly Logical. Should the assembly summary file be retrieved and stored?
#'   #' If true, the assembly summaries are retrieved and stored at the assembly_path. If false,
#'   #' the assembly summary at the assembly_path is read and used.
#'   #' @param file type at the assembly_path. Either 'tsv' or 'csv'
#'   if(updateAssembly == T)
#'   {
#'     # Over 700,000
#'     assembly_k_genbank <- getKingdomAssemblySummary("genbank")
#'
#'     ## Over 200,000
#'     assembly_k_refseq <- getKingdomAssemblySummary("refseq")
#'
#'     assembly_all <- bind_rows(assembly_k_genbank,assembly_k_refseq) #%>%
#'                       #select(all_of(keep))
#'
#'     switch(file_type,
#'            "tsv" = write_tsv(assembly_all, assembly_path),
#'            "csv" = write_csv(assembly_all,assembly_path)
#'     )
#'
#'   }
#'   else
#'   {
#'     assembly_all <- switch(file_type,
#'            "tsv" =  read_tsv(assembly_path),
#'            "csv" =  read_csv(assembly_path)
#'            )
#'   }
#'
#'   keep = c("assembly_accession", "taxid", "species_taxid", "organism_name")
#'
#'   assembly_all <- assembly_all %>% select(all_of(keep))
#'
#'
#' }
