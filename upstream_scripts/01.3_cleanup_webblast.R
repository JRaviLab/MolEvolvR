# Script to cleanup the raw web BLAST output
# adds column names; calculates PcPositive w.r.t. query

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# compute cvm location for acc2lin.R
source("/data/research/jravilab/molevol_scripts/R/lineage.R")
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")

## Read in data file path as a string
args <- commandArgs(trailingOnly = TRUE)

## clean up function
# takes in path to blast result file as a string
cleanup_webblast <- function(infile_blast) {

  # column names for blast results
  web_blast_cols <- c("Query", "AccNum",
                      "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                      "QStart", "QEnd", "SStart", "SEnd",
                      "EValue", "BitScore", "PcPosOrig")

  # read in blast results from input variable, set column names to blast_cols (created above)
  blast_out <- read_tsv(file = infile_blast, col_names = web_blast_cols, skip = 7)

  cleanedup_blast <- blast_out %>%
    # remove extra characters/names from sseqid, sscinames, and staxids columns
    mutate(AccNum = gsub('\\|', '', AccNum)) %>%
    mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
    mutate(PcIdentity = round(PcIdentity, 2)) # %>%
    # normalize percent positive by multiplying the original ppos by the length of the subject protein
    #   length and dividing by the query protein length
    # mutate(PcPositive = round(x = (PcPosOrig * AlnLength/QLength), digits = 2))
  
  wblast_lins <- add_lins(df = cleanedup_blast, assembly_path = '/data/research/jravilab/common_data/assembly_summary_genbank.txt', lineagelookup_path = '/data/research/jravilab/common_data/lineage_lookup.txt') %>%
     add_name()
  
  # create name for new output file to be created
  file_name <- gsub(pattern = '.txt', replacement = '', x = infile_blast) %>%
    paste0('.cln.txt')
  # write the cleaned up data to new file
  write_tsv(cleanedup_blast, file_name, col_names = T)
}

cleanup_webblast(args[1])

