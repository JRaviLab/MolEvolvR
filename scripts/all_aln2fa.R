## To generate fasta files from aln files
## Pre-requisites to generate MSA and Phylogenetic Tree
## List of functions: generate_all_aln2fa | part of pre-msa-tree.R

## Modified: Dec 24, 2019
## Authors: Janani Ravi (@jananiravi)

#################
## Pkgs needed ##
#################
library(here); library(tidyverse)
library(data.table)
#library(seqRFLP)
conflicted::conflict_prefer("filter", "dplyr")

##################
## The Function ##
##################
generate_all_aln2fa <- function(aln_path=here("data/rawdata_aln/"),
                                fa_outpath=here("data/alns/"),
                                lin_file=here("data/rawdata_tsv/all_semiclean.txt"),
                                reduced=F)
{
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to all alignment files w/ accessions & DAs?
  #' @param aln_path Character. Path to alignment files.
  #' Default is 'here("data/rawdata_aln/")'
  #' @param lin_file Character. Path to file. Master protein file with AccNum & lineages.
  #' Default is 'here("data/rawdata_tsv/all_semiclean.txt")'
  #' @param fa_outpath Character. Path to the written fasta file.
  #' Default is 'here("data/alns/")'.
  #' @param reduced Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
  #' Default is 'FALSE'.
  #' @examples generate_all_aln2fa()
  #' @details The alignment files would need two columns separated by spaces: 1. AccNum and 2. alignment. The protein homolog file should have AccNum, Species, Lineages.
  #' @note Please refer to the source code if you have alternate + file formats and/or column names.

  library(here)
  library(tidyverse)
  # aln_path <- here("data/rawdata_aln/")
  # outpath <- here("data/alns/")
  # lin_file <- here("data/rawdata_tsv/all_semiclean.txt")

  aln_filenames <- list.files(path=aln_path, pattern="*.aln")
  aln_filepaths <- paste0(aln_path, aln_filenames)
  variable <- str_replace_all(basename(aln_filenames),
                              pattern=".aln", replacement="")

  ## Using purrr::pmap
  aln2fa_args <- list(aln_file=aln_filepaths,
                      fa_outpath=paste0(fa_outpath, variable, ".fa"))
  pmap(.l=aln2fa_args, .f=convert_aln2fa,
       lin_file=lin_file,
       reduced=reduced)
}

### UNUSED ###
## Using lapply | has issues
# lapply(X=aln_filepaths,
#        FUN=convert_aln2fa,
#        aln_file=aln_filepaths,
#        lin_file=rep(lin_file, length(aln_files)),
#        fa_outpath=paste0(fa_outpath, variable, ".fa"),
#        reduced=rep(reduced, length(aln_files)))

## Using a simple for loop
# for(i in 1:length(aln_filepaths)){
#   print(i)
#   convert_aln2fa(aln_file=aln_filepaths[i],
#                  fa_outpath=paste0(fa_outpath, variable[i], ".fa"),
#                  lin_file=lin_file,
#                  reduced=reduced)
# }


## Individual file format
# convert_aln2fa(aln_file = "data/rawdata_aln/pspa_snf7.gismo.aln",
#                 lin_file = lin_file,
#                reduced = reduced, fa_outpath = "data/alns/pspa_snf7.fa")
