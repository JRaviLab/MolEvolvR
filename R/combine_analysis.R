# Script to combine full analysis and IPR analysis files
# Uses R/combine_files.R and writes out cln_combined.tsv and ipr_combined.tsv
# Created: Dec 19, 2020
# Lauren Sosinski, Janani Ravi

source('R/combine_files.R')
source(here('R/combine_files.R')) # to work locally

combine_full <- function(inpath, ret=FALSE) {
  ## Combining full_analysis files
  full_combnd <- combine_files(inpath,
                               pattern="*.full_analysis.tsv", skip=0,
                               col_names=T)

  write_tsv(x=full_combnd, col_names=T,
            file=paste0(inpath, "/", "cln_combined.tsv"))
  if(ret)
  {
    return(full_combnd)
  }

}

combine_ipr <- function(inpath, ret=FALSE) {
  ## Combining clean ipr files
  ipr_combnd <- combine_files(inpath,
                              pattern="*.iprscan_cln.tsv", skip=0,
                              col_names=T)

  write_tsv(x=ipr_combnd, col_names=T,
            file=paste0(inpath, "/", "ipr_combined.tsv"))
  if(ret)
  {
    return(ipr_combnd)
  }
}
