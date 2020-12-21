library(tidyverse)
library(data.table)
library(furrr)

# source lineage.R
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")
source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")

## load files in
# ipr2da function
ipr2da <- function(infile_ipr, acc2info, prefix, analysis=c("Pfam","SMART", "CDD", "TIGRFAM",
                                                  "Phobius", "Gene3D", "TMHMM", "SignalP_EUK",
                                                  "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE"))
{
 
#infile_ipr <- 'sample5.iprscan.tsv'
#acc2info <- 'sample.acc2info.tsv'
#prefix <- 'sample5'

 # read in lookup table
  lookup_tbl <- fread("/data/research/jravilab/common_data/cln_lookup_tbl.tsv", header = T, fill = T) %>%
    arrange() %>% distinct()

  # read in iprscan results
  ipr_in <- read_tsv(infile_ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))

  # read in acc2info
  acc2info_out <- read_tsv(file = acc2info, col_names = T) %>%
    mutate(FullAccNum = gsub('\\|', '', FullAccNum)) %>%
    mutate(AccNum = gsub('.*[a-z]', '', FullAccNum)) %>%
    select(-AccNum.noV, -FullAccNum)

  prot_in_da <- ipr_in %>%
    merge(lookup_tbl, by = "DB.ID", all.x = T)
  
  # split dataframe into unique proteins
  x <- split(x = prot_in_da, f = prot_in_da$AccNum)

  plan(strategy = "multicore", .skip = T)

  # within each data.table
  domarch <- future_map(x, function(y) {
    acc_row <- data.frame(AccNum = y$AccNum[1],  stringsAsFactors = F)
    DAs <- data.frame(matrix(nrow = 1, ncol = length(analysis) ))
    DA <- y %>% group_by(Analysis) %>% arrange(StartLoc)
    i = 1
    for(a in analysis) {
      a_da <- DA %>% filter(Analysis == a)
      if (a == "SignalP_EUK" || a == "SignalP_GRAM_NEGATIVE" || a == "SignalP_GRAM_POSITIVE") {
        var_shortname = "DB.ID" }
      else {
        var_shortname = "ShortName" }
      var_shortname_sym = sym(var_shortname)
      a_da <- a_da %>%
        ungroup() %>%
        select({{var_shortname_sym}}) %>%
        filter(!is.na({{var_shortname_sym}})) %>%
        filter(!is.null({{var_shortname_sym}})) %>%
        pull(var_shortname) %>%
        paste(collapse = "+")
      DAs[1,i] = a_da
      i=(i+1)
    }

    colnames(DAs) = paste("DomArch", analysis, sep = ".")
    return(cbind(acc_row, DAs))
  })
  
  # read in lineage mapping file
  lineage_map <- fread("/data/research/jravilab/common_data/lineage_lookup.txt", header = T, fill = T)
 
  # combine domarchs to one data frame, merge w/ acc2info
  domarch2 <- do.call(rbind.data.frame, domarch) %>%
    merge(acc2info_out, by = "AccNum")
  
  # add lineage to domarch, remove extra species column
  domarch_lins <- merge(domarch2, lineage_map, by = "TaxID", all.x = T) %>%
   select(-Species.x)
  
  # change Species.y colname to Species
  names(domarch_lins)[names(domarch_lins) == 'Species.y'] <- 'Species'
  # add name column to domarch+lineage dataframe
  domarch_lins <- domarch_lins %>% add_name()

  # save domarch_lins file
  write_tsv(domarch_lins, file = paste0(prefix, '.ipr_domarch.tsv'), append = F, na = 'NA')
  
  # return domarch2 dataframe to append to blast results if given
  return(domarch2)
}

## function for adding results from ipr2da to blast results

append_ipr <- function(ipr_da, blast, prefix) {
  #ipr_domarch <- read_tsv(ipr_da, col_names = T)
  blast_out <- read_tsv(blast, col_names = T)
  
  blast_ipr <- merge(blast_out, ipr_da, by = 'AccNum')
  write_tsv(blast_ipr, file = paste0(prefix, '.full_analysis.tsv'), na = 'NA')
}

# accept command line args
args <- commandArgs(trailingOnly=TRUE)

## perform ipr2da on iprscan results
da <- ipr2da(infile_ipr = args[1], acc2info =  args[2], prefix = args[3])

## if blast results are provided, call append_ipr
if (is.null(args[4]) | is.na(args[4])) {
   print('No blast results provided, moving on.') 
} else {
   append_ipr(ipr_da=da, blast=args[4], prefix=args[3]) 
}


