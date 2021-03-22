# Created by Lauren Sosinski
# Last edit on 2021.03.22
# Script used to create the lookup table used by MolEvolvR in ipr2da

library(tidyverse)
library(data.table)

# read in data tables w/o column names
cdd_cols <- c('ID', 'DB.ID', 'ShortName', 'Description', 'Misc.')
cdd<- fread('cddid_all.tbl', sep = '\t', col.names = cdd_cols, quote = '\"') %>%
  select(-Misc.) %>%
  select(DB.ID, ShortName, Description, ID)
cdd$ID <- as.character(cdd$ID)

pfam_cols <- c('DB.ID', 'ID', 'Misc', 'ShortName', 'Description')
pfam <- read_tsv('Pfam-A.clans.tsv', col_names = pfam_cols) %>%
  select(-Misc)
pfam$ID <- as.character(pfam$ID)

# read in others
tmhmm <- read_tsv('TMHMM_lookup.tsv')
sigp <- read_tsv('SignalP_lookup.tsv')
phobius <- read_tsv('Phobius_lookup.tsv')
panther <- read_tsv('pnthr_db_cln.tsv')
g3d <- read_tsv('ln-cath.tsv')

# merge data frames
lookup_cln <- cdd %>%
  bind_rows(tmhmm, sigp, phobius, panther, pfam, g3d)

# write the file
write_tsv(lookup_cln, "cln_lookup_tbl.tsv", col_names = T)
