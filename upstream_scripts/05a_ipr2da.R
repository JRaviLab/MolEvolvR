library(tidyverse)
library(data.table)
library(furrr)

# source lineage.R
source("R/pre-msa-tree.R")
source("R/colnames_molevol.R")

# ipr2da function
ipr2da <- function(infile_ipr, prefix,
                   analysis=c("Pfam","SMART","Phobius",
                              "Gene3D", "TMHMM", "PANTHER", "ProSiteProfiles",
                              "SUPERFAMILY", "MobiDBLite", "TIGRFAM"))
{
  # read in cleaned up iprscan results
  ipr_in <- read_tsv(infile_ipr, col_names = T)

  # split dataframe into unique proteins
  x <- split(x = ipr_in, f = ipr_in$AccNum)

  plan(strategy = "multicore", .skip = T)

  # within each data.table
  #domarch <- map(x, function(y) {
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

  # select relevant rows from ipr input to add to domarch
  ipr_select <- ipr_in %>%
    select(Name, AccNum, Species, TaxID, Lineage, ProteinName, SourceDB, Completeness) %>%
    distinct()

  # combine domarchs to one data frame, merge w/ acc2info
  domarch2 <- do.call(rbind.data.frame, domarch)

  domarch_lins <- domarch2 %>%
    merge(ipr_select, by = 'AccNum', all.x = T)

  # save domarch_lins file
  write_tsv(domarch_lins, file = paste0(prefix, '.ipr_domarch.tsv'),
            append = F, na = 'NA')

  # return domarch2 dataframe to append to blast results if given
  return(domarch2)
}

## function to add results from ipr2da to blast results
append_ipr <- function(ipr_da, blast, prefix) {
  blast_out <- read_tsv(blast, col_names = T)

  blast_ipr <- merge(blast_out, ipr_da, by = 'AccNum', all.x = T) %>%
    mutate(Species = Species.y, TaxID = TaxID.y, Name = Name.y, Lineage = Lineage.x) %>%
    select(-Species.x, -Species.y, -TaxID.x, -TaxID.y, -Name.x, -Name.y, -Lineage.x, -Lineage.y)

  write_tsv(blast_ipr, file = paste0(prefix, '.full_analysis.tsv'), na = 'NA')
}

# accept command line args
args <- commandArgs(trailingOnly=TRUE)

## perform ipr2da on iprscan results
da <- ipr2da(infile_ipr = args[1], prefix = args[2])

## if blast results are provided, call append_ipr
if (is.null(args[3]) | is.na(args[3])) {
   print('No blast results provided, moving on.')
} else {
   append_ipr(ipr_da=da, blast=args[3], prefix=args[2])
}
