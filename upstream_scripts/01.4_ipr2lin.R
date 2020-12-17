library(tidyverse)
library(data.table)

# source lineage script for add_name and add_lins
source('/data/research/jravilab/molevol_scripts/R/lineage.R')
source('/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R')
source('/data/research/jravilab/molevol_scripts/R/colnames_molevol.R')

# add lineage to iprscan results

ipr <- '../laurensosinski/data/WP_001901328_full/WP_001901328.iprscan.tsv'
acc2tax <- '../laurensosinski/data/WP_001901328_full/WP_001901328.acc2tax.tsv'
suffix <- 'WP_001901328'

ipr2lin <- function(ipr, acc2tax, suffix) {
  # read in iprscan results,
  ipr_in <- read_tsv(ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))

  a2t_out <- fread(input = acc2tax, sep = '\t', header = T) %>%
    mutate(FullAccNum = gsub('\\|', '', FullAccNum)) %>%
    mutate(FullAccNum = gsub('.*[a-z]', '', FullAccNum))

  # run add_lins function on ipr output
  ipr_tax <- merge(ipr_in, a2t_out, by.x = "AccNum", by.y = "FullAccNum")
  lineage_map <- fread("/data/research/jravilab/common_data/lineage_lookup.txt", header = T, fill = T)
  # # get lineage path as argument, it'll be changed depending on who is running it
  # # have default argument also for where the files are located
  ipr_lin <- merge(ipr_tax, lineage_map, by = "TaxID", all.x = T)

  # write results to file
  write_tsv(ipr_lin, file = paste0(suffix, '.iprscan_lins.tsv'))
}

## load files in
args <- commandArgs(trailingOnly = TRUE)

## call function
ipr2lin(args[1], args[2])
