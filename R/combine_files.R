## Simple function to combine multiple files of the same type
## Args: file path, search pattern, column names
## Created: Nov 23, 2020

## !!YET TO ROXYGENIZE!!

##################
## Libs & Paths ##
##################
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(here))

# STARTING DIRECTORY
# molevol_scripts
# source("R/colnames_molevol.R")
# source(here("R/colnames_molevol.R")) # to work locally
# inpath <- c("../molevol_data/project_data/slps/full_analysis_20201207/")

#################################
## Fn to combine similar files ##
#################################
#' Download the combined assembly summaries of genbank and refseq
#'
#' @author Janani Ravi
#'
#' @param inpath String of 'master' path where the files reside (recursive=T)
#' @param pattern Character vector containing search pattern for files
#' @param delim
#' @param skip
#' @param col_names Takes logical T/F arguments OR column names vector;
#' usage similar to col_names parameter in `readr::read_delim`
#'
#' @importFrom purrr pmap_dfr
#' @importFrom readr cols
#'
#' @return
#' @export
#'
#' @examples
combine_files <- function(inpath = c("../molevol_data/project_data/phage_defense/"),
    pattern = "*full_analysis.tsv",
    delim = "\t", skip = 0,
    col_names = T) {
    source_files <- dir(
        path = inpath, pattern = pattern,
        recursive = T
    )
    source_files_path <- paste0(inpath, source_files)
    combnd_files <- source_files_path %>%
        list() %>%
        pmap_dfr(read_delim,
            delim = delim,
            col_names = col_names,
            col_types = cols(Score = col_character()), # to avoid datatype error
            skip = skip,
            .id = "ByFile"
        )

    return(combnd_files)
}

# pmap_dfr(fread, fill=T, na.strings=c(""), header=T) %>%
# rbindlist(., use.names=T, fill=T, idcol=NULL)

#################
## Sample Runs ##
#################
# ## Combining full_analysis files
# full_combnd <- combine_files(inpath,
#                             pattern="*full_analysis.txt", skip=0,
#                             col_names=T)
#
# write_tsv(x=full_combnd, col_names=T,
#           path="../molevol_data/project_data/slps/full_combined.tsv")
#
# ## Combining clean files
# cln_combnd <- combine_files(inpath,
#                             pattern="^.*cln.txt", skip=0,
#                             col_names=T)
#
# write_tsv(x=cln_combnd, col_names=T,
#           path="../molevol_data/project_data/slps/cln_combined.tsv")
#
#
# ## Less helpful examples!
# ## Combining BLAST files
# ## Likely makes no sense since clustering is done per query
# cl_blast_combnd <- combine_files(inpath,
#                                  pattern="^.*refseq.1e-5.txt", skip=0,
#                                  col_names=cl_blast_colnames) %>%
#   select(-PcPositive, -ClusterID)
#
# ## Combining IPR files
# ## Likely makes no sense since there may be repeated AccNum from indiv. files!
# ipr_combnd <- combine_files(inpath,
#                             pattern="*iprscan.lins*",  skip=0,
#                             col_names=ipr_colnames)
#
# write_tsv(x=ipr_combnd, col_names=T,
#           path="../molevol_data/project_data/slps/ipr_combined.tsv")
