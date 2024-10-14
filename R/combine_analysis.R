# Script to combine full analysis and IPR analysis files
# Uses R/combine_files.R and writes out cln_combined.tsv and ipr_combined.tsv
# Created: Dec 19, 2020
# Lauren Sosinski, Janani Ravi

# source("R/combine_files.R")
# source(here("R/combine_files.R")) # to work locally

#' Combining full_analysis files
#'
#' @param inpath
#' @param ret
#'
#' @importFrom readr write_tsv
#'
#' @return
#' @export
#'
#' @examples
combineFullAnalysis <- function(inpath, ret = FALSE) {
    ## Combining full_analysis files
    full_combnd <- combineFiles(inpath,
        pattern = "*.full_analysis.tsv", skip = 0,
        col_names = T
    )

    write_tsv(
        x = full_combnd, col_names = T,
        file = paste0(inpath, "/", "cln_combined.tsv")
    )
    if (ret) {
        return(full_combnd)
    }
}

#' Combining clean ipr files
#'
#' @param inpath
#' @param ret
#'
#' @importFrom readr write_tsv
#'
#' @return
#' @export
#'
#' @examples
combineIPR <- function(inpath, ret = FALSE) {
    ## Combining clean ipr files
    ipr_combnd <- combineFiles(inpath,
        pattern = "*.iprscan_cln.tsv", skip = 0,
        col_names = T
    )

    write_tsv(
        x = ipr_combnd, col_names = T,
        file = paste0(inpath, "/", "ipr_combined.tsv")
    )
    if (ret) {
        return(ipr_combnd)
    }
}
