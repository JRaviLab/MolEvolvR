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
#' @param inpath Character. The master directory path where the files reside.
#' The search is recursive (i.e., it will look in subdirectories as well).
#' @param pattern Character. A search pattern to identify files to be combined.
#' Default is "*full_analysis.tsv".
#' @param delim Character. The delimiter used in the input files.
#' Default is tab ("\t").
#' @param skip Integer. The number of lines to skip at the beginning of each file.
#' Default is 0.
#' @param col_names Logical or character vector. If TRUE, the first row of each file
#' is treated as column names. Alternatively, a character vector can
#' be provided to specify custom column names.
#'
#' @importFrom purrr pmap_dfr
#' @importFrom readr cols col_character
#'
#' @return A data frame containing the combined contents of all matched files.
#' Each row will include a new column "ByFile" indicating the source file of the data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' combined_data <- combineFiles(inpath = "../molevol_data/project_data/phage_defense/")
#' }
combineFiles <- function(inpath = c("../molevol_data/project_data/phage_defense/"),
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
# full_combnd <- combineFiles(inpath,
#                             pattern="*full_analysis.txt", skip=0,
#                             col_names=T)
#
# write_tsv(x=full_combnd, col_names=T,
#           path="../molevol_data/project_data/slps/full_combined.tsv")
#
# ## Combining clean files
# cln_combnd <- combineFiles(inpath,
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
# cl_blast_combnd <- combineFiles(inpath,
#                                  pattern="^.*refseq.1e-5.txt", skip=0,
#                                  col_names=cl_blast_colnames) %>%
#   select(-PcPositive, -ClusterID)
#
# ## Combining IPR files
# ## Likely makes no sense since there may be repeated AccNum from indiv. files!
# ipr_combnd <- combineFiles(inpath,
#                             pattern="*iprscan.lins*",  skip=0,
#                             col_names=ipr_colnames)
#
# write_tsv(x=ipr_combnd, col_names=T,
#           path="../molevol_data/project_data/slps/ipr_combined.tsv")
