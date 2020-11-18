library(tidyverse)
library(data.table)
library(furrr)

# compute cvm location for acc2lin.R
#source("/data/research/jravilab/laurensosinski/scripts/acc2lin.R")
# local location for acc2lin.R
source("scripts/acc2lin.R")

## load file in
# on HPCC
args <- commandArgs(trailingOnly = TRUE)

# locally
# args <- c("data/molevolvr_outputs/WP_020839904.1_Vibrio_parahaemolyticus_out/WP_020839904.1_Vibrio_parahaemolyticus.iprscan.tsv",
#           "data/molevolvr_outputs/WP_020839904.1_Vibrio_parahaemolyticus_out/WP_020839904.1_Vibrio_parahaemolyticus.refseq.1e-5.cln.txt")

# ipr2da function
ipr2da <- function(infile_ipr, infile_blast, suffix, analysis=c("Pfam","SMART", "CDD", "TIGRFAM",
                                                        "Phobius", "Gene3D", "TMHMM", "SignalP_EUK",
                                                        "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE"))
  {
  infile_ipr = "data/molevolvr_outputs/slps/ADY20190.1_EA1_Bacillus_thuringiensis_out/ADY20190.1_EA1_Bacillus_thuringiensis.iprscan.tsv"
  infile_blast = "data/molevolvr_outputs/slps/ADY20190.1_EA1_Bacillus_thuringiensis_out/ADY20190.1_EA1_Bacillus_thuringiensis.refseq.1e-5.cln.txt"
  analysis=c("Pfam","SMART", "CDD", "TIGRFAM", "Phobius", "Gene3D", "TMHMM", "SignalP_EUK",
             "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE")
  suffix = "ADY20190.1_EA1_Bacillus_thuringiensis"

  # creating column names for input files and lookup tables
  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate", "IPRAcc", "IPRDesc")

  lookup_cath <- fread('/data/research/jravilab/common_data/cln-cath.tsv') %>%
    arrange() %>% distinct()

  lookup_tbl <- fread("/data/research/jravilab/common_data/lookup_tbl.tsv", header = T, fill = T) %>%
    arrange() %>% distinct() %>%
    bind_rows(lookup_cath)

  # read in iprscan and blast results
  blast_out <- fread(infile_blast, header = T, keepLeadingZeros = T)

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
        var_shortname = "Short_Name" }
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

  # TaxID to lineage
  blast_out$TaxID <- as.integer(blast_out$TaxID)
  lineage_map <- fread("data/ReferenceFiles/lineagelookup.txt", sep = "\t")
  # get lineage path as argument, it'll be changed depending on who is running it
  # have default argument also for where shit is
  mergedLins <- merge(blast_out, lineage_map, by.x = "TaxID", by.y="tax_id", all.x = T)
  updated_blast <- merge(mergedLins, domarch2, by = "AccNum")

  #write_tsv(mergedLins, file = paste0(suffix, '.iprscan_lins.tsv'))
  write_tsv(updated_blast, file = infile_blast)

  # # add DA col + localization prediction col
  ## pfam column, COG column, TIGR, Superfam, SMART, localization
  ## COGS in RPSBLAST
}

ipr2da(args[1],args[2])
