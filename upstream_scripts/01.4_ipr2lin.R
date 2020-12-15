library(tidyverse)
library(data.table)

# source lineage script for add_name and add_lins
source('/data/research/jravilab/molevol_scripts/R/lineage.R')
source('/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R')
source('/data/research/jravilab/molevol_scripts/R/colnames_molevol.R')

# add lineage to iprscan results

ipr2lin <- function(ipr, suffix) {
  # ipr column names
  #ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
  #                  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
  #                  "Status", "RunDate", "IPRAcc", "IPRDesc")

  # read in iprscan results,
  ipr_in <- read_tsv(ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))
  
  # run add_lins function on ipr output
  ipr_lin <- ipr_in %>%
    add_lins(assembly_path = '/data/research/jravilab/common_data/assembly_summary_refseq.txt',
             lineagelookup_path = '/data/research/jravilab/common_data/lineage_lookup.txt') %>%
    add_name()
  
  # write results to file
  write_tsv(ipr_lin, file = paste0(suffix, '.iprscan_lins.tsv'))
}

## load files in
args <- commandArgs(trailingOnly = TRUE)

## call function 
ipr2lin(args[1], args[2])
