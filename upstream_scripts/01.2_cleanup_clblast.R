# Script to cleanup the raw BLAST script
# adds column names; calculates PcPositive w.r.t. query

library(tidyverse)
library(data.table)

# compute cvm location for acc2lin.R
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")

## Read in data file path as a string
args <- commandArgs(trailingOnly = TRUE)

## clean up function
# takes in path to blast result file as a string
cleanup_clblast <- function(infile_blast) {

  # column names for blast results
  cl_blast_cols <- c("Query", "SAccNum", "AccNum", "SAllSeqID", "STitle", "Species", "TaxID",
                     "PcIdentity", "AlnLength", "Mismatch", "Gapopen","QStart", "QEnd",
                     "QLength", "SStart", "SEnd", "SLength", "EValue", "BitScore",
                     "PcPosOrig")

  # read in blast results from input variable, set column names to blast_cols (created above)
  blast_out <- read_tsv(file = infile_blast, col_names = cl_blast_cols)

  cleanedup_blast <- blast_out %>%
    # remove extra characters/names from sseqid, sscinames, and staxids columns
    mutate(AccNum = gsub('\\|', '', AccNum)) %>%
    mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
    mutate(Species = gsub(';.*$', '', Species)) %>%
    mutate(TaxID = gsub(';.*$', '', TaxID)) %>%
    mutate(PcIdentity = round(PcIdentity, 2)) %>%
    # normalize percent positive by multiplying the original ppos by the length of the subject protein
    #   length and dividing by the query protein length
    mutate(PcPositive = round(x = (PcPosOrig * AlnLength/QLength), digits = 2))

  # # TaxID to lineage
  cleanedup_blast$TaxID <- as.integer(cleanedup_blast$TaxID)
  lineage_map <- fread("/data/research/jravilab/common_data/lineagelookup.txt", sep = "\t")
  # # get lineage path as argument, it'll be changed depending on who is running it
  # # have default argument also for where the files are located
  mergedLins <- merge(cleanedup_blast, lineage_map, by.x = "TaxID", by.y="TaxID", all.x = T)
  
  blast_names <- add_name(mergedLins)

  # create name for new output file to be created
  file_name <- gsub(pattern = '.txt', replacement = '', x = infile_blast) %>%
    paste0('.cln.txt')
  # write the cleaned up data to new file
  write_tsv(blast_names, file_name, col_names = T)
}

cleanup_clblast(args[1])
