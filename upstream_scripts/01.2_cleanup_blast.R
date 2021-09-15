# Script to cleanup the raw BLAST script
# adds column names; calculates PcPositive w.r.t. query

library(tidyverse)
library(data.table)

# source scripts w/ add_name() and colnames
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")
source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")

## Read in data file path as a string
args <- commandArgs(trailingOnly=TRUE)


## clean up function
# takes in path to blast result file as a string
cleanup_blast <- function(infile_blast, acc2info, prefix, wblast) {

  # read in acc2info file, cleanup FullAccNum column
  acc2info_out <- fread(input = acc2info, sep = "\t", header = T, fill = T) %>%
    mutate(FullAccNum = gsub('\\|', '', FullAccNum)) %>%
    mutate(FullAccNum = gsub('.*[a-z]', '', FullAccNum))

  # read in blast results, set colnames, cleanup results, merge acc2info output
  if (wblast == "T") {
     blast_out <- fread(input = infile_blast, sep = '\t', header = F,
       col.names = web_blastp_hit_colnames, fill = T)
            query <- blast_out[1,]$Query
     cleanedup_blast <- blast_out %>%
     # remove extra characters/names from sseqid, sscinames, and staxids columns
       mutate(AccNum = gsub('\\|', '', AccNum)) %>%
       mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
       mutate(PcIdentity = round(as.double(PcIdentity), 2))
     # merge blast out with acc2info out
     cleanedup_blast <- merge(cleanedup_blast, acc2info_out, by.x = "AccNum", by.y = "FullAccNum")
     # find query in acc2info, extract & set Length as QLength
     query_prot <- as.character(cleanedup_blast[1,2])
     q_row <- acc2info_out[ grepl(query_prot, FullAccNum) ]
     qlen <- as.character(q_row[1,Length])
     # if row in acc2info contains query_prot, extract Length value & populate column with it
     cleanedup_blast$QLength <- qlen
     cleanedup_blast <- cleanedup_blast %>%
       mutate(PcPosOrig = as.numeric(PcPosOrig)) %>%
       mutate(AlnLength = as.numeric(AlnLength)) %>%
       mutate(QLength = as.numeric(QLength)) %>%
       mutate(PcPositive = round(x = (PcPosOrig * AlnLength/QLength), digits = 2))

  } else if (wblast == "F") {
     blast_out <- read_tsv(file = infile_blast, col_names = cl_blast_colnames)
     query <- blast_out[1,]$Query
     cleanedup_blast <- blast_out %>%
       # remove extra characters/names from sseqid, sscinames, and staxids columns
       mutate(AccNum = gsub('\\|', '', AccNum)) %>%
       mutate(AccNum = gsub('.*[a-z]', '', AccNum)) %>%
       mutate(Species = gsub(';.*$', '', Species)) %>%
       mutate(PcIdentity = round(PcIdentity, 2)) %>%
       # normalize percent positive by multiplying the original ppos by the length of the subject protein
       #   length and dividing by the query protein length
       mutate(PcPositive = round(x = (PcPosOrig * AlnLength/QLength), digits = 2))

     # merge blast out with acc2info out
     cleanedup_blast <- merge(cleanedup_blast, acc2info_out, by.x = "AccNum", by.y = "FullAccNum", all = TRUE) %>%
       select(-Species.x, -TaxID.x)
     names(cleanedup_blast)[names(cleanedup_blast) == 'Species.y'] <- 'Species'
  }

  # # TaxID to lineage
  cleanedup_blast$TaxID <- as.integer(cleanedup_blast$TaxID)
  lineage_map <- fread("/data/research/jravilab/common_data/lineage_lookup.txt", header = T, fill = T)
  # # get lineage path as argument, it'll be changed depending on who is running it
  # # have default argument also for where the files are located
  mergedLins <- merge(cleanedup_blast, lineage_map, by = "TaxID", all.x = T) %>%
    mutate(Species=Species.y, Spp.blast=Species.x) %>%
    select(any_of(cl_blast_postcln_cols))

  blast_names <- add_name(mergedLins)
  print(blast_names)
  # begin query name addition
  query_row <- subset(blast_names, AccNum==prefix)
  query_name <- query_row$Name
  if(!is.null(query_name)){
    query_name <- query
  }
  blast_names$QueryName <- query_name
  blast_names <- subset(blast_names, Query!="NA")

  # create name for new output file to be created
  file_name <- paste0(prefix, '.blast.cln.tsv')
  # write the cleaned up data to new file
  write_tsv(blast_names, file_name, col_names = T)
}

#cleanup_clblast(infile_blast, wblast)
cleanup_blast(args[1], args[2], args[3], args[4])
