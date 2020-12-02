library(tidyverse)
library(data.table)
library(furrr)

# source acc2lin.R
source("/data/research/jravilab/molevol_scripts/upstream_scripts/acc2lin.R")

## load files in
args <- commandArgs(trailingOnly = TRUE)

# ipr2da function
ipr2da <- function(infile_ipr, infile_blast, analysis=c("Pfam","SMART", "CDD", "TIGRFAM",
                                                        "Phobius", "Gene3D", "TMHMM", "SignalP_EUK",
                                                        "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE")) #, suffix)
  {
  # creating column names for input files and lookup tables
  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate", "IPRAcc", "IPRDesc")

  lookup_tbl <- fread("/data/research/jravilab/common_data/cln-lookup_tbl.tsv", header = T, fill = T) %>%
    arrange() %>% distinct()

  # read in blast results
  blast_out <- fread(infile_blast, header = T, keepLeadingZeros = T)

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
        pull(var_shortname) %>%
        paste(collapse = "+")
      DAs[1,i] = a_da
      i=(i+1)
    }

    colnames(DAs) = paste("DomArch", analysis, sep = ".")
    return(cbind(acc_row, DAs))
  })

  domarch2 <- do.call(rbind.data.frame, domarch)

  ## extract accession numbers, sort by unqiue accs
  # accs <- data.frame(ipr_in$AccNum) %>% unique()
  ## get lineage using Sam's acc2lin function
  # acc2lin(accs, "../common_data/assembly_summary2020-11-02.txt",
  #         "../common_data/lineagelookup.txt")
  ## create new cleanedup iprscan file
  # write_tsv(mergedLins, file = paste0(suffix, '.iprscan_lins.tsv'))

  # TaxID to lineage
  blast_out$TaxID <- as.integer(blast_out$TaxID)
  lineage_map <- fread("/data/research/jravilab/common_data/lineagelookup.txt", sep = "\t")
  # get lineage path as argument, it'll be changed depending on who is running it
  # have default argument also for where things are
  mergedLins <- merge(blast_out, lineage_map, by.x = "TaxID", by.y="tax_id", all.x = T)
  updated_blast <- merge(mergedLins, domarch2, by = "AccNum")

  lin <- select(blast, AccNum, TaxID, Species, Lineage)

  ipr_lin <- merge(ipr_in, mergedLins, by = 'AccNum')

  write_tsv(ipr_lin, paste0(suffix, 'iprscan.lins.tsv', collapse = '.'), append = FALSE)
  #write_tsv(mergedLins, file = paste0(suffix, '.iprscan_lins.tsv'))
  write_tsv(updated_blast, file = infile_blast, append = F)
}

ipr2da(args[1],args[2])
