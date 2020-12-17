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

infile_blast = "../laurensosinski/data/WP_001901328_full/WP_001901328.refseq.1e-5.txt"
acc2tax = "../laurensosinski/data/WP_001901328_full/WP_001901328.acc2tax.tsv"
wblast = "T"


## clean up function
# takes in path to blast result file as a string
cleanup_blast <- function(infile_blast, acc2tax, wblast) {

  # read in acc2tax file, cleanup FullAccNum column
  a2t_out <- fread(input = acc2tax, sep = '\t', header = T) %>%
    mutate(FullAccNum = gsub('\\|', '', FullAccNum)) %>%
    mutate(FullAccNum = gsub('.*[a-z]', '', FullAccNum))

  # read in blast results, set colnames, cleanup results, merge acc2taxlin output
  if (wblast == "T") {
     blast_out <- read_tsv(file = infile_blast, col_names = web_blastp_desc_colnames)
     cleanedup_blast <- blast_out %>%
       # remove extra characters/names from sseqid, sscinames, and staxids columns
       mutate(AccNum = gsub('\\|', '', AccNum)) %>%
       mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
       mutate(PcIdentity = round(PcIdentity, 2))

     # merge blast out with acc2tax out
     cleanedup_blast <- merge(cleanedup_blast, a2t_out, by.x = "AccNum", by.y = "FullAccNum") %>%
       select(-Species.x, -TaxID.x)
     names(cleanedup_blast)[names(cleanedup_blast) == 'Species.y'] <- 'Species'
     names(cleanedup_blast)[names(cleanedup_blast) == 'TaxID.y'] <- 'TaxID'
  } else if (wblast == "F") {
     blast_out <- read_tsv(file = infile_blast, col_names = cl_blast_colnames)
     cleanedup_blast <- blast_out %>%
       # remove extra characters/names from sseqid, sscinames, and staxids columns
       mutate(AccNum = gsub('\\|', '', AccNum)) %>%
       mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
       mutate(Species = gsub(';.*$', '', Species)) %>%
       mutate(PcIdentity = round(PcIdentity, 2)) %>%
       # normalize percent positive by multiplying the original ppos by the length of the subject protein
       #   length and dividing by the query protein length
       mutate(PcPositive = round(x = (PcPosOrig * AlnLength/QLength), digits = 2))

     # merge blast out with acc2tax out
     cleanedup_blast <- merge(cleanedup_blast, a2t_out, by.x = "AccNum", by.y = "FullAccNum") %>%
       select(-Species.x, -TaxID.x)
     names(cleanedup_blast)[names(cleanedup_blast) == 'Species.y'] <- 'Species'
     names(cleanedup_blast)[names(cleanedup_blast) == 'TaxID.y'] <- 'TaxID'
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
cleanup_blast(args[1], args[2], args[3])
