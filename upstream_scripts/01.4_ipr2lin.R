library(tidyverse)
library(data.table)

# add lineage to iprscan results

ipr2lin <- function(ipr, blast, suffix) {
  # ipr column names
  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate", "IPRAcc", "IPRDesc")

  # read in blast results
  blast_out <- fread(infile_blast, header = T, keepLeadingZeros = T)

  # read in iprscan results,
  ipr_in <- read_tsv(infile_ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))
  
  # extract necessary columns from blast output file to add to ipr results
  lin <- select(blast_out, AccNum, TaxID, Species, Lineage)
  
  # combine ipr file w/ lineage columns from blast
  # run add_name function to add species names to ipr results via taxID
  ipr_lin <- merge(ipr_in, lin, by = 'AccNum', all.x = T) %>%
    add_name()
  
  # write results to file
  write_tsv(ipr_lin, file = paste0(suffix, '.iprscan_lins.tsv'))
}

## load files in
args <- commandArgs(trailingOnly = TRUE)

## call function 
ipr2lin(args[1], args[2], args3[])
