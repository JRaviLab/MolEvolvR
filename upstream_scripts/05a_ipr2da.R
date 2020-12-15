library(tidyverse)
library(data.table)
library(furrr)

# source lineage.R
# source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")

## load files in
args <- commandArgs(trailingOnly = TRUE)

# ipr2da function
ipr2da <- function(infile_ipr, suffix, analysis=c("Pfam","SMART", "CDD", "TIGRFAM",
                                                  "Phobius", "Gene3D", "TMHMM", "SignalP_EUK",
                                                  "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE"))
  {
  # creating column names for input files and lookup tables
  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate", "IPRAcc", "IPRDesc")

  lookup_tbl <- fread("/data/research/jravilab/common_data/cln_lookup_tbl.tsv", header = T, fill = T) %>%
    arrange() %>% distinct()

  # read in blast results
  #blast_out <- fread(infile_blast, header = T, keepLeadingZeros = T)
  
  # read in iprscan results,
  ipr_in <- read_tsv(infile_ipr, col_names = ipr_colnames) %>%
    mutate(DB.ID = gsub('G3DSA:', '', DB.ID))

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
        var_shortname = "Short.Name" }
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

  domarch2 <- do.call(rbind.data.frame, domarch)
  
  write_tsv(domarch2, file = paste0(suffix, '.ipr_domarch.tsv'), append = F, na = 'NA')
  
  return(domarch2)
}

## function for adding results from ipr2da to blast results

append_ipr <- function(ipr_da, blast, suffix) {
  ipr_domarch <- fread(ipr_da, header = T, fill = T)
  blast_out <- fread(blast, header = T, fill = T)
  
  blast_ipr <- merge (blast_out, ipr_domarch, by = 'AccNum')
  write_tsv(blast_ipr, file = paste0(suffix, '.cln.clust.ipr.tsv'), na = 'NA')
}

## perform ipr2da on iprscan results
da <- ipr2da(args[1],args[2])

## if blast results are provided, call append_ipr
if (is.null(args[3]) | length(args[3] == 0)) {
  print('No blast results provided, moving on.') 
} else if (args[3] == grep('.cln.clust.txt', args[3])) {
  append_ipr(da, blast=args[3], suffix=args[2]) 
}


