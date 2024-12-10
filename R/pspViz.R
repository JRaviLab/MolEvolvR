# figure 4, from the PSP paper, https://journals.asm.org/doi/epub/10.1128/msystems.00847-23
# code from https://github.com/JRaviLab/psp_app/blob/main/psp_process.Rmd

library(tidyverse)
library(data.table)

# source("/data/research/jravilab/molevol_scripts/R/cleanup.R")  # in package
# source("/data/research/jravilab/molevol_scripts/R/summarize.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/plotting.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/reverse_operons.R") # in package
# source("/data/research/jravilab/molevol_scripts/R/networks_domarch.R") # in package

cln_colnames <- c("AccNum", "DomArch", "GenContext", "Length",  "GeneName",
                  "Lineage", "Species", "GCA_ID", "GeneDesc",
                  "DomArch.repeats", "GenContext.repeats")

query_domains <- read_delim("data/acc_files/query_domains.txt",
                            delim="\t", col_names=T)
lineage_rename <- read_tsv("data/acc_files/lineage_rename.txt",
                           col_names=T)

data_cleanup <- function(all_raw) {
  # Cleanup Species
  all_cln <- all_raw %>% cleanup_species(remove_empty=FALSE)
  # Cleanup Gene Description
  all_cln <- all_cln %>% cleanup_GeneDesc(column="GeneDesc") # Remove trailing '.'
  # CLeanup Lineage
  all_cln <- all_cln %>% cleanup_lineage(lins_rename=lineage_rename)

  # Cleanup DomArch
  # For the repeats
  all_cln <- all_cln %>%
    cleanup_domarch(domains_rename=domains_rename,
                    repeat2s=FALSE)
                    # remove_tails=F,      # FIXME # keep if it exists & works for DomArch
                    # remove_empty=F)      # FIXME # keep if it exists & works for DomArch


  all_cln$DomArch.repeats <- all_cln$DomArch

  # Removing duplicate AccNum w/ different DomArchs
  all_cln$DomArch.uncompressed <- all_cln$DomArch.repeats
  # !! repeat2s: deprecation notice for funs & list
  all_cln <- repeat2s(all_cln, "DomArch.uncompressed")

  # Extract unique rows
  all_cln <- all_cln %>% distinct()
  # Pick longer of the duplicated AccNum DomArchs
  all_cln <- all_cln %>% 
    pick_longer_duplicate("DomArch.uncompressed")

  all_cln <- all_cln %>%
    cleanup_domarch(domains_rename=domains_rename,
                    domains_keep=NULL,   
                    domains_ignore=NULL, 
                    repeat2s=TRUE,
                    remove_tails=F,      
                    remove_empty=F)      

  # FIXME # Not using genomic context parts of the script, yet
  # Cleanup GenContext
  # Calls reverseOperons
  #  all_cln <- all_cln %>%
  #    cleanup_gencontext(domains_rename=domains_rename,
  #                      repeat2s=FALSE)
  #
  #  all_cln$GenContext.repeats <- all_cln$GenContext
  #
  #  all_cln <- all_cln %>%
  #    cleanup_gencontext(domains_rename=domains_rename,
  #                      repeat2s=T)

  # Subset essential columns for downstream analyses
  all_cln <- all_cln %>% select(all_of(cln_colnames))

  # Remove repeated observations
  all_cln <- all_cln %>% distinct()




  prot <- all_cln
}

plotLineageBarStack <- function(job_dir) {
  # FIXME # figure out if all_raw can be read from the job folder
  all_w_extra <- read_tsv(file="data/rawdata_tsv/all_with_extrapspasnf7.tsv")
  all_raw <- all_w_extra

  # cleanup function on all_raw to produce prot
  # FIXME # Check if this is needed
  prot <- data_cleanup(all_raw)

  # rename lineages
  prot$Lineage.reduced = prot$Lineage %>%
   str_replace(pattern = "^eukaryota.*", replacement = "eukaryota") %>%
   str_replace(pattern = "^viruses.*", replacement = "viruses")

  prot$Lineage = prot$Lineage.reduced

  # Include only DAs ≥ min.cutoff
  prot_cutoff <- prot %>% total_counts(cutoff=100)

  prot_cutoff <- prot_cutoff %>% filter(totalcount>=50)

  cutoff_perc <- 100 - prot_cutoff$CumulativePercent[nrow(prot_cutoff)]

  stacked_lin <- stacked_lin_plot(prot= prot, column="DomArch",
                                cutoff=cutoff_perc,
                                label.size=20, 
                                legend.position=c(0.6, 0.40), legend.text.size=20,
                                legend.size=1)

  return(stacked_lin)

}
