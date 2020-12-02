library(tidyverse)
library(data.table)

# read in files
args <- commandArgs(trailingOnly = TRUE)

rps2da <- function(infile_rps, infile_blast) {
  # analyses to filter by
  analysis <- c("COG") #, "PRK", "cd")
  # column names for rpsblast output
  rps_colnames <- c("AccNum", "DB.ID", "DBSeqID",
                    "PcIdentity.Dom", "PcPosOrig.Dom", #"PcPos.Dom",
                    "AlnLength", "Mismatch",
                    "SStart", "SEnd", "DStart", "DEnd",
                    "EValue", "BitScore") #, "TaxID") # TaxID missing (NA); remove?
  # reading in rpsblast output
  rpsout <- read_tsv(file = infile_rps, col_names = rps_colnames)
  # reading in blast output
  blast_out <- read_tsv(file = infile_blast, col_names = T)
  # lookup table
  lookup <- read_tsv("/data/research/jravilab/common_data/lookup_tbl.tsv")
  # clean up rpsblast cdd id number
  rpsout <- rpsout %>%
    mutate(DB.ID = gsub(pattern = "CDD:", replacement = "", x=DB.ID)) %>%
    mutate(DB.ID = as.character(DB.ID))

  rps_cdd <- merge(rpsout, lookup, by.x="DB.ID", by.y ="ID")

  #### filter by analysis ####
  # filter by COG
  rps_cog <- rps_cdd %>% filter(grepl(pattern = "COG", x = rps_cdd$DB.ID.y))

  # split data frame by accession numbers
  x <- split(x = rps_cog, f = rps_cog$AccNum)

  domarch <- map(x, function(y) {
    acc_row <- data.frame(AccNum = y$AccNum[1],  stringsAsFactors = F)
    DAs <- data.frame(matrix(nrow = 1, ncol = 1))
    DA <- y %>% arrange(SStart)
    i = 1

    a_da <- DA %>%
      select(Short.Name) %>%
      filter(!is.na(Short.Name)) %>%
      pull(Short.Name) %>%
      paste(collapse = "+")
    DAs[1,i] = a_da
    i=(i+1)

    colnames(DAs) = paste("DomArch", "COG", sep = ".")
    return(cbind(acc_row, DAs))
  })

  domarch2 <- do.call(rbind.data.frame, domarch)
  rps_da <- merge(blast_out, domarch2)

  #### SAVE RPS W/ CDD DATA TABLE ####
  write_tsv(rps_da, file = infile_blast)

}

rps2da(args[1], args[2])
