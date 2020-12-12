library(tidyverse)
library(rentrez)
library(future)
library(furrr)
library(data.table)


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

  assembly_all <- assembly_all %>% data.table::setnames(
    old = c("taxid","refseq_category","species_taxid","organism_name","infraspecific_name","genome_rep"),
    new = c( "TaxID", "RefseqCategory","Parent.TaxID","Species","Spp.Strain","GenomeStatus"),
    skip_absent = T)

  #
  # dplyr::rename("AssemblyID"="assembly_accession",
  #                                              "TaxID" = "taxid",
  #                                              "RefseqCategory" = "refseq_category",
  #                                              "Parent.TaxID" = "species_taxid",
  #                                              "Species" = "organism_name",
  #                                              "Spp.Strain" = "infraspecific_name",
  #                                              "GenomeStatus" = "genome_rep")

  fwrite(assembly_all, outpath, sep = "\t")
}




# Go from the GCA_ID column to tax IDs using the assembly file
GCA2lin <- function(prot_data, assembly_path = "data/acc_files/assembly_summary20201018.txt",
                     lineagelookup_path = "data/lineagelookup.txt", acc_col = "Protein" )
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
  assembly_summary <- setnames(assembly_summary, "AssemblyID", "GCA_ID")

  mergedTax <- merge(x = prot_data,y = assembly_summary,by = "GCA_ID", all.x = T)

  accessions = prot_data %>% pull(acc_col) %>% unique()

  best_rows = integer(length(accessions))
  for(i in 1:length(accessions))
  {
    # browser()
    acc = accessions[i]
    acc_inds = which(mergedTax$Protein == acc)
    if(length(acc_inds) == 1)
    {
      best_rows[i] = acc_inds
    }
    else if(length(acc_inds) > 1)
    {
      # refseq inds take precedence
      refseq_inds = acc_inds[which(mergedTax[acc_inds,]$Source == "RefSeq")]
      if(length(refseq_inds) == 1)
      {
        # only one refseq ind, keep that
        best_rows[i] = refseq_inds[1]
        next
      }
      else if( length(refseq_inds) > 1)
      {
        # more than 1 refseq ind
        refseq_complete = refseq_inds[which(mergedTax[refseq_inds,]$assembly_level == "Complete Genome")]
        if(length(refseq_complete) >= 1)
        {
          # Take the first one with refseq complete
          best_rows[i] = refseq_complete[1]
          next
        }
        # take the first refseq if no complete
        best_rows[i] = refseq_inds[1]
      }
      else
      {
        # only non refseq: should I just choose the first one?
        if(any(mergedTax[acc_inds]$assembly_level == "Complete Genome"))
        {
          best_rows[i] = acc_inds[which(mergedTax[acc_inds]$assembly_level == "Complete Genome")][1]
        }
        else{
          best_rows[i] = acc_inds[1]
        }
      }
    }
  }

  mergedTax = mergedTax[best_rows,]

  lineage_map <- fread(lineagelookup_path, sep = "\t")
  lineage_map = lineage_map[,!"Species"]

  mergedLins <- merge(mergedTax, lineage_map, by.x = "TaxID", by.y="TaxID",
                      all.x = T)

  return(mergedLins)
}


# https://stackoverflow.com/questions/18730491/sink-does-not-release-file
sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}


add_lins <- function(df, acc_col = "AccNum", assembly_path,
                     lineagelookup_path, ipgout_path = NULL)
{
  s_acc_col = sym(acc_col)
  accessions = df %>% pull(acc_col)
  lins = acc2lin(accessions, assembly_path, lineagelookup_path, ipgout_path)

  # Drop a lot of the unimportant columns for now? will make merging much easier
  lins <- lins[,c("Strand","Start","Stop", "Nucleotide Accession", "Source",
                  "Id", "Strain"):= NULL]
  lins <- unique(lins)

  # dup <- lins %>% group_by(Protein) %>% summarize(count = n()) %>% filter(count > 1) %>%
  #   pull(Protein)

  merged = merge(df, lins, by.x = acc_col, by.y = "Protein", all.x = TRUE)
  return(merged)
}


acc2lin <- function(accessions,  assembly_path, lineagelookup_path,ipgout_path = NULL )
{
  #'@author Samuel Chen, Janani Ravi
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

  # if(tmp_ipg)
  # {
  #   unlink(tempdir(), recursive = T)
  # }

  # cols = c("TaxID","GCA_ID", "Protein", "Protein Name", "Species", "Lineage")
  # lins = unique(lins[,..cols])

  return(lins)
}



efetch_ipg <- function(accNum_vec, out_path)
{
  #'@author Samuel Chen, Janani Ravi
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
      if(x%%9 == 0)
      {
        Sys.sleep(1)
      }
      cat(
        # entrez_fetch(id = partitioned_acc[[x]],
        #              db = "protein",
        #              rettype = "txt",# parsed = T,
        #              api_key = "55120df9f5dddbec857bbb247164f86a2e09"
        entrez_fetch(id = partitioned_acc[[x]],
                     db = "ipg",
                     rettype = "xml",# parsed = T,
                     api_key = "55120df9f5dddbec857bbb247164f86a2e09"
        )
      )
    })
    sink(NULL)
  }
}

ipg2lin <- function(accessions, ipg_file, assembly_path, lineagelookup_path)
{
  #'@author Samuel Chen, Janani Ravi
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

  lins <- GCA2lin(prot_data = ipg_dt, assembly_path, lineagelookup_path)
  lins <- lins[!is.na(Lineage)] %>% unique()

  return(lins)
}


