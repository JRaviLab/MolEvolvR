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

domains_rename <- read_delim("data/acc_files/domains_rename.txt",
                             delim="\t", col_names=TRUE)

# domains_ignore <- read_delim("data/acc_files/domains_ignore.txt",
#                              delim="\t", col_names=T)

domains_keep <- read_delim("data/acc_files/domains_keep.txt",
                           delim="\t", col_names=T)

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
                    domains_keep=NULL,   # filter applied to only ClustName for now.
                    domains_ignore=NULL, #!! should check and remove soon!
                    repeat2s=FALSE,
                    remove_tails=F,      #new! check below if it works!
                    remove_empty=F)      #new! needed?

  all_cln$DomArch.repeats <- all_cln$DomArch

  # Removing duplicate AccNum w/ different DomArchs
  all_cln$DomArch.uncompressed <- all_cln$DomArch.repeats
  # !! repeat2s: deprecation notice for funs & list
  all_cln <- repeat2s(all_cln, "DomArch.uncompressed",
                      excluded_prots=c("PspC", "LiaI-LiaF-TM"))

  # Extract unique rows
  all_cln <- all_cln %>% distinct()
  # Pick longer of the duplicated AccNum DomArchs
  all_cln <- all_cln %>% pick_longer_duplicate("DomArch.uncompressed")

  all_cln <- all_cln %>%
    cleanup_domarch(domains_rename=domains_rename,
                    domains_keep=NULL,   
                    domains_ignore=NULL, 
                    repeat2s=TRUE,
                    remove_tails=F,      
                    remove_empty=F)      

  # Cleanup GenContext
  # Calls reverse_operons
  all_cln <- all_cln %>%
    cleanup_gencontext(domains_rename=domains_rename,
                      repeat2s=FALSE)

  all_cln$GenContext.repeats <- all_cln$GenContext

  all_cln <- all_cln %>%
    cleanup_gencontext(domains_rename=domains_rename,
                      repeat2s=T)

  # Subset essential columns for downstream analyses
  all_cln <- all_cln %>% select(all_of(cln_colnames))

  # Remove repeated observations
  all_cln <- all_cln %>% distinct()



  # Replace SIG+Toastrack with TM+Toastrack
  all_cln$DomArch <- all_cln$DomArch %>% str_replace_all(pattern = "SIG\\+Toastrack", replacement = "TM+Toastrack")
  all_cln$DomArch.repeats <- all_cln$DomArch.repeats %>% str_replace_all(pattern = "SIG\\+Toastrack", replacement = "TM+Toastrack")

  all_cln$GenContext <- all_cln$GenContext %>% str_replace_all(pattern = "SIG\\+Toastrack", replacement = "TM+Toastrack")
  all_cln$GenContext.repeats <- all_cln$GenContext.repeats %>% str_replace_all(pattern = "SIG\\+Toastrack", replacement = "TM+Toastrack")

  # Replace SIG+DUF4178 with TM+4178
  all_cln$DomArch <- all_cln$DomArch %>% str_replace_all(pattern = "SIG\\+DUF4178", replacement = "TM+DUF4178")
  all_cln$DomArch.repeats <- all_cln$DomArch.repeats %>% str_replace_all(pattern = "SIG\\+DUF4178", replacement = "TM+DUF4178")

  all_cln$GenContext <- all_cln$GenContext %>% str_replace_all(pattern = "SIG\\+DUF4178", replacement = "TM+DUF4178")
  all_cln$GenContext.repeats <- all_cln$GenContext.repeats %>% str_replace_all(pattern = "SIG\\+DUF4178", replacement = "TM+DUF4178")

  # Replace SIG+PspA with TM+PspA
  all_cln$DomArch <- all_cln$DomArch %>% str_replace_all(pattern = "SIG\\+PspA", replacement = "TM+PspA")
  all_cln$DomArch.repeats <- all_cln$DomArch.repeats %>% str_replace_all(pattern = "SIG\\+PspA", replacement = "TM+PspA")

  all_cln$GenContext <- all_cln$GenContext %>% str_replace_all(pattern = "SIG\\+PspA", replacement = "TM+PspA")
  all_cln$GenContext.repeats <- all_cln$GenContext.repeats %>% str_replace_all(pattern = "SIG\\+PspA", replacement = "TM+PspA")


  # Convert Sig+Snf7 to TM+Snf7
  all_cln$DomArch <- all_cln$DomArch %>% str_replace_all(pattern = "SIG\\+Snf7", replacement = "TM+Snf7")
  all_cln$DomArch.repeats <- all_cln$DomArch.repeats %>% str_replace_all(pattern = "SIG\\+Snf7", replacement = "TM+Snf7")
  all_cln$GenContext <- all_cln$GenContext %>% str_replace_all(pattern = "SIG\\+Snf7", replacement = "TM+Snf7")
  all_cln$GenContext.repeats <- all_cln$GenContext.repeats %>% str_replace_all(pattern = "SIG\\+Snf7", replacement = "TM+Snf7")


  ## Write cleaned up file 
  # write_tsv(all_cln, "data/rawdata_tsv/all_clean_combined_20210329.txt")


  prot <- all_cln
}

psp_fig4 <- function(job_dir) {
  # FIXME: figure out if all_raw can be read from the job folder
  all_w_extra <- read_tsv(file="data/rawdata_tsv/all_with_extrapspasnf7.tsv")
  # all_w_extra$GI
  all_raw <- all_w_extra

  # similar to the Rmd, we run the cleanup function on all_raw to produce prot
  prot <- data_cleanup(all_raw)

  # the rest of the figure generation follows as-is
  prot$Lineage.reduced = prot$Lineage %>% str_replace(pattern = "^eukaryota.*", replacement = "eukaryota") %>% str_replace(pattern = "^viruses.*", replacement = "viruses")

  prot$Lineage = prot$Lineage.reduced

  nopspa <- prot %>%
    filter_by_doms(column="DomArch", 
                  doms_remove=c("PspA", "PspA(s)", "Snf7"),
                  ignore.case=T)

  # Include only DAs ≥ min.cutoff
  nopspa_tc <- nopspa %>% total_counts(cutoff=100)

  nopspa_tc <- nopspa_tc %>% filter(totalcount>=50)

  cutoff_perc <- 100 - nopspa_tc$CumulativePercent[nrow(nopspa_tc)]

  stlin_nops <- stacked_lin_plot(prot=nopspa, column="DomArch",
                                cutoff=cutoff_perc,
                                label.size=20, 
                                legend.position=c(0.6, 0.40), legend.text.size=20,
                                legend.size=1)

  # res =  stacked_lin_plot(prot=nopspa, column="DomArch",
  #                                cutoff=cutoff_perc,
  #                                label.size=20, 
  #                                legend.position=c(0.6, 0.40), legend.text.size=20,
  #                                legend.size=1)

  return(stlin_nops)

  # FA: writing out the file commented out in favor of returning the object

  # Decrease width, increase height, increase font size
  # 
  # ggsave("stackedLinPlot_50tc.png", stlin_nops,
  #       width=15,
  #       height=18,
  #       dpi=300)
}
