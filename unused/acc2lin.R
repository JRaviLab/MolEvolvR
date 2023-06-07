suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(biomartr))

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

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
# https://stackoverflow.com/questions/18730491/sink-does-not-release-file
sink.reset <- function(){
for(i in seq_len(sink.number())){
    sink(NULL)
  }
}
acc2lin <- function(accessions,  assembly_path, lineagelookup_path,ipgout_path = NULL )
{
#'@author Samuel Chen
#'@description This function combines 'efetch_ipg()' and 'ipg2lin()' to map a set
#'of protein accessions to their assembly (GCA_ID), tax ID, and lineage.
#'@param accessions Character vector of protein accessions
#'@param assembly_path String of the path to the assembly_summary path
#'This file can be generated using the "DownloadAssemblySummary()" function
#'@param lineagelookup_path String of the path to the lineage lookup file
#'(taxid to lineage mapping). This file can be generated using the
#'@param ipgout_path Path to write the results of the efetch run of the accessions
#'on the ipg database. If NULL, the file will not be written. Defaults to NULL
tmp_ipg = F
if(is.null(ipgout_path))
  {
tmp_ipg = T
ipgout_path = tempfile("ipg", fileext =".txt")
  }
  efetch_ipg(accessions, out_path= ipgout_path )
lins <- ipg2lin(accessions, ipgout_path, assembly_path, lineagelookup_path)
if(tmp_ipg)
  {
    unlink(tempdir(), recursive = T)
  }
return(lins)
}
efetch_ipg <- function(accNum_vec, out_path)
{
#'@author Samuel Chen
#'@description Perform efetch on the ipg database and write the results to out_path
#'@param accNum_vec Character vector containing the accession numbers to query on
#'the ipg database
#'@param out_path Path to write the efetch results to
if(length(accNum_vec) > 0){
partition <- function(v, groups){
# Partition data to limit number of queries per second for rentrez fetch:
# limit of 10/second w/ key
l <- length(v)
partitioned <- list()
for(i in 1:groups)
      {
partitioned[[i]] <- v[seq.int(i,l,groups)]
      }
return(partitioned)
    }
    plan(strategy = "multiprocess", .skip = T)
min_groups = length(accNum_vec)/200
groups <- min(max(min_groups,15) ,length(accNum_vec))
partitioned_acc <- partition(accNum_vec, groups )
    sink(out_path)
a <- future_map(1:length(partitioned_acc), function(x)
    {
# Avoid hitting the rate API limit
if(x%%10 == 0)
      {
        Sys.sleep(1)
      }
        cat(
          entrez_fetch(id = partitioned_acc[[x]],
db = "ipg",
rettype = "xml",
api_key = "YOUR_KEY_HERE"
          )
        )
    })
    sink(NULL)
  }
}
ipg2lin <- function(accessions, ipg_file, assembly_path, lineagelookup_path)
{
#'@author Samuel Chen
#'@description Takes the resulting file of an efetch run on the ipg database and
#'append lineage, and taxid columns
#'@param accessions Character vector of protein accessions
#'@param ipg_file Filepath to the file containing results of an efetch run on the
#'ipg database. The protein accession in 'accessions' should be contained in this
#'file
#'@param assembly_path String of the path to the assembly_summary path
#'This file can be generated using the "DownloadAssemblySummary()" function
#'@param lineagelookup_path String of the path to the lineage lookup file
#'(taxid to lineage mapping). This file can be generated using the
#'"create_lineage_lookup()" function
ipg_dt <- fread(ipg_file, sep = "\t", fill = T)
ipg_dt <- ipg_dt[Protein %in% accessions]
ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")
lins <- GCA2Lins(prot_data = ipg_dt, assembly_path, lineagelookup_path)
return(lins)
}

suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(here))
# suppressPackageStartupMessages(library(biomartr))


create_lineage_lookup <- function(lineage_file = "data/rankedlineage.dmp",
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
  #'Choices include: "supperkingdom", "phylum",   "class","order", "family",
  #'"genus", and "species"

  shorten_NA <- function(Lineage)
  {
    first_NA = str_locate(Lineage, "NA")[1]
    if(is.na(first_NA) )
    {
      # No NAs
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
    mutate(Lineage = unlist(map(Lineage, shorten_NA)))



  write_tsv(x = rankedLinsCombined, path = outfile)
}
