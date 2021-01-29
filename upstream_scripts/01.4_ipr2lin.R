library(tidyverse)
library(data.table)

# source lineage script for add_name and add_lins
# /data/research/jravilab/molevol_scripts/
source('/data/research/jravilab/molevol_scripts/R/lineage.R')
source('/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R')
source('/data/research/jravilab/molevol_scripts/R/colnames_molevol.R')

# add lineage to iprscan results

ipr2lin <- function(ipr, acc2info, prefix) {

  # read in iprscan results
  # duplicate rows in iprscan file
  ipr_in <- read_tsv(ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))

  acc2info_out <- fread(input = acc2info, sep = '\t', header = T, fill = T) %>%
    mutate(FullAccNum = gsub('\\|', '', FullAccNum)) %>%
    mutate(FullAccNum = gsub('.*[a-z]', '', FullAccNum))

  # merge ipr file with acc2info file
  ipr_in <- ipr_in %>%
    mutate(AccNum.noV = gsub('\\.[0-9]', '', AccNum))

  ipr_tax <- merge(ipr_in, acc2info_out, by = 'AccNum.noV', all.x = T)

  # read in lineage map
  lineage_map <- fread("/data/research/jravilab/common_data/lineage_lookup.txt", header = T, fill = T)
  #lineage_map <- fread("../ReferenceFiles/lineage_lookup.txt", header = T, fill = T)

  # merge ipr+info w/ lineage, remove extra species column
  ipr_lin <- merge(ipr_tax, lineage_map, by = "TaxID", all.x = T) %>%
    mutate(Species = Species.y) %>%
    select(-Species.x, -Species.y)

  # add lookup table to iprscan file
  lookup_tbl <- fread(input = '/data/research/jravilab/common_data/cln_lookup_tbl.tsv', sep = '\t', header = T, fill = T) %>%
  #lookup_tbl <- fread(input = '../ReferenceFiles/cln_lookup_tbl.tsv', sep = '\t', header = T, fill = T) %>%
    distinct()

  # run add_name f(x) on ipr+lineage dataframe
  ipr_lin <- ipr_lin %>% add_name() %>%
    mutate(Name = gsub('^_', '', Name))

  # add domarch info to iprscan + lineage df, only keep what's in x
  ipr_cln <- merge(ipr_lin, lookup_tbl, by = 'DB.ID', all.x = T, all.y = F)

  # populate empty description/short name columns
  for (i in 1:nrow(ipr_cln)) {
    if ((is.na(ipr_cln$ShortName[i]) || ipr_cln$ShortName[i] == '') && (is.na(ipr_cln$SignDesc[i]) || ipr_cln$SignDesc[i] == '-')) {
      ipr_cln$SignDesc[i] = ipr_cln$IPRDesc[i]
      if (length(ipr_cln$LookupTblDesc[i]) != 0) {
        ipr_cln$ShortName[i] = ipr_cln$LookupTblDesc[i]
      }
    }
    if (is.na(ipr_cln$ShortName[i]) || ipr_cln$ShortName[i] == '') {
      ipr_cln$ShortName[i] = ipr_cln$SignDesc[i]
    }
    if (is.na(ipr_cln$SignDesc[i]) || ipr_cln$SignDesc[i] == '-') {
      ipr_cln$SignDesc[i] = ipr_cln$ShortName[i]
    }
  }

  # rename unclear/duplicated columns
  names(ipr_cln)[names(ipr_cln) == 'Description.x'] <- 'ProteinName'
  names(ipr_cln)[names(ipr_cln) == 'Description.y'] <- 'LookupTblDesc'

  # create label column to use in ipr2viz
  ipr_cln <- ipr_cln %>%
    mutate(Label = strtrim(ShortName, 30)) %>%
    mutate(Label = gsub(', .*', '', Label)) %>%
    mutate(Label = gsub('C-terminal region of a signal', 'C-term signal peptide', Label)) %>%
    mutate(Label = gsub('N-terminal region of a signal', 'N-term signal peptide', Label)) %>%
    mutate(Label = gsub('Twin arginine translocation \\(T', 'Tat signal profile', Label)) %>%
    mutate(Label = gsub('GLUTATHIONE HYDROLASE PROENZYM', 'GLUTATHIONE HYDROLASE PROENZYME', Label)) %>%
    mutate(Label = gsub('N-terminal nucleophile aminohy', 'Ntn hydrolases', Label)) %>%
    mutate(Label = gsub('Region of a membrane-bound pro', 'cytoplasmic reg of mem-bound prot', Label)) %>%
    mutate(ShortName = gsub('Region of a membrane-bound protein predicted to be outside the membrane,
                            in the cytoplasm.', 'cytoplasmic reg of mem-bound prot', ShortName))

  # write results to file
  write_tsv(ipr_cln, file = paste0(prefix, '.iprscan_cln.tsv'))
}

## load files in
args <- commandArgs(trailingOnly = TRUE)

## call function
ipr2lin(args[1], args[2], args[3])
