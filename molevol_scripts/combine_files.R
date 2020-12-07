## Simple function to combine multiple files of the same type
## Args: file path, search pattern, column names
## Created: Nov 23, 2020

## !!YET TO ROXYGENIZE!!

##################
## Libs & Paths ##
##################
library(tidyverse)
library(data.table)
library(here)

# STARTING DIRECTORY
# molevol_scripts
source(here("molevol_scripts/colnames_molevol.R"))
inpath <- c("../molevol_data/project_data/slps/full_analysis_20201202/")

#################################
## Fn to combine similar files ##
#################################
combine_files <- function(inpath=c("../molevol_data/project_data/phage_defense/"),
                          pattern="^WP_.*refseq.1e-5.txt",
                          col_names=cl_blast_colnames)
{
  #' Download the combined assembly summaries of genbank and refseq
  #' @author Janani Ravi
  #' @param inpath String of 'master' path where the files reside (recursive=T)
  #' @param pattern Character vector containing search pattern for files of a kind
  #' @param col_names Takes logical T/F arguments OR column names vector; usage similar to col_names parameter in `readr::read_delim`

  source_files <- dir(path=inpath, pattern=pattern,
                      recursive=T)
  source_files_path <- paste0(inpath, source_files)
  combnd_files <- source_files_path %>% list %>%
    # pmap_dfr(fread, fill=T, na.strings=c(""), header=T)
    pmap_dfr(read_tsv, col_names=col_names, .id="ByFile")
  return(combnd_files)
}

#################
## Sample Runs ##
#################
## Combining clean files
cln_combnd <- combine_files(inpath,
                            pattern="^WP_.*cln.*",
                            col_names=T)
# cln_combnd_noda <- cln_combnd %>%
#   select(-starts_with("DomArch"), -tax_name, -Lineage)

write_tsv(x=cln_combnd, col_names=T,
          path="../molevol_data/project_data/slps/cln_combined.tsv")


## Less helpful examples!
## Combining BLAST files
## Likely makes no sense since clustering is done per query
cl_blast_combnd <- combine_files(inpath,
                                 pattern="^WP_.*refseq.1e-5.txt",
                                 col_names=cl_blast_colnames) %>%
  select(-PcPositive, -ClusterID)

## Combining IPR files
## Likely makes no sense since there may be repeated AccNum from indiv. files!
ipr_combnd <- combine_files(inpath,
                            pattern="*iprscan.lins*",
                            col_names=ipr_colnames)

write_tsv(x=ipr_combnd, col_names=T,
          path="../molevol_data/project_data/slps/ipr_combined.tsv")
