# Script to cleanup the raw BLAST script
# adds column names; calculates PcPositive w.r.t. query

library(tidyverse)
library(data.table)

# compute cvm location for acc2lin.R
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")
#source("R/pre-msa-tree.R")
source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")
#source("R/colnames_molevol.R")

## Read in data file path as a string
args <- commandArgs(trailingOnly = TRUE)

#infile_blast = "../laurensosinski/data/molevolvr_outputs/saureus_runs/sausa300_0200_out/sausa300_0200.refseq.1e-5.txt"
#wblast = T


## clean up function
# takes in path to blast result file as a string
cleanup_clblast <- function(infile_blast, acc2tax, wblast) {

  # read in blast results from input variable, set column names to blast_cols (created above)
  if (wblast == "T") {
     blast_out <- read_tsv(file = infile_blast, col_names = web_blastp_desc_colnames)
     cleanedup_blast <- blast_out %>%
       # remove extra characters/names from sseqid, sscinames, and staxids columns
       mutate(AccNum = gsub('\\|', '', AccNum)) %>%
       mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
       mutate(TaxID = gsub(';.*$', '', TaxID)) %>%
       mutate(PcIdentity = round(PcIdentity, 2))
  } else if (wblast == "F") {
     blast_out <- read_tsv(file = infile_blast, col_names = cl_blast_colnames)
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
  }

  # # TaxID to lineage
  cleanedup_blast$TaxID <- as.integer(cleanedup_blast$TaxID)
  lineage_map <- fread("/data/research/jravilab/common_data/lineage_lookup.txt", header = T, fill = T)
  # # get lineage path as argument, it'll be changed depending on who is running it
  # # have default argument also for where the files are located
  mergedLins <- merge(cleanedup_blast, lineage_map, by = "TaxID", all.x = T) %>%
    mutate(Species=Species.y, Spp.blast=Species.x) %>%
    select(all_of(cl_blast_postcln_cols))

  blast_names <- add_name(mergedLins)

  # create name for new output file to be created
  file_name <- gsub(pattern = '.txt', replacement = '', x = infile_blast) %>%
    paste0('.cln.tsv')
  # write the cleaned up data to new file
  write_tsv(blast_names, file_name, col_names = T)
}

#cleanup_clblast(infile_blast, wblast)
cleanup_clblast(args[1], args[2], args[3])
