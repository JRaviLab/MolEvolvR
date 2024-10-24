# Script to combine full analysis and IPR analysis files
# Uses R/combine_files.R and writes out cln_combined.tsv and ipr_combined.tsv
# Created: Dec 19, 2020
# Lauren Sosinski, Janani Ravi

# source("R/combine_files.R")
# source(here("R/combine_files.R")) # to work locally

#' Combining full_analysis files
#'
#' @param inpath Character. The path to the directory containing the 
#' `.full_analysis.tsv` files to be combined.
#' @param ret Logical. If TRUE, the function will return the combined data frame. 
#' Default is FALSE, meaning it will only write the file and not return the data.
#'
#' @importFrom readr write_tsv
#'
#' @return If `ret` is TRUE, a data frame containing the combined data from all 
#' input files. If `ret` is FALSE, the function writes the combined data to a 
#' TSV file named `cln_combined.tsv` in the specified directory and returns NULL.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' combined_data <- combineFullAnalysis("path/to/full_analysis/files", ret = TRUE)
#' }
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
#' @param inpath Character. The path to the directory containing the 
#' `.iprscan_cln.tsv` files to be combined.
#' @param ret Logical. If TRUE, the function will return the combined data frame. 
#' Default is FALSE, meaning it will only write the file and not return the data.
#'
#' @importFrom readr write_tsv
#'
#' @return If `ret` is TRUE, a data frame containing the combined data from all 
#' input files. If `ret` is FALSE, the function writes the combined data to a 
#' TSV file named `ipr_combined.tsv` in the specified directory and returns NULL.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' combineIPR <- combine_ipr("path/to/ipr/files", ret = TRUE)
#' }
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
