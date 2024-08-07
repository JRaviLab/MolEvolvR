# suppressPackageStartupMessages(library(tidyverse))

#' Run DELTABLAST to find homologs for proteins of interest
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @param deltablast_path
#' @param db_search_path Path to the BLAST databases
#' @param db
#' @param query
#' @param evalue
#' @param out
#' @param num_alignments
#' @param num_threads
#'
#' @return
#' @export
#'
#' @examples
run_deltablast <- function(deltablast_path, db_search_path,
                           db = "refseq", query, evalue = "1e-5",
                           out, num_alignments, num_threads = 1) {
  #' Run DELTABLAST to find homologs for proteins of interest
  #' @author Samuel Chen, Janani Ravi
  #' @param path
  #' @param db_search_path Path to the BLAST databases
  #' @param db
  #' @param query
  #' @param evalue
  #' @param out
  #' @param num_threads
  #' @param num_alignments

  start <- Sys.time()

    system(paste0("export BLASTDB=/", db_search_path))

    system2(
        command = deltablast_path,
        args = c(
            "-db", db,
            "-query", query,
            "-evalue", evalue,
            "-out", out,
            "-num_threads", num_threads,
            "-num_alignments", num_alignments
            #   ,"-outfmt", outfmt
        )
    )
    print(Sys.time() - start)
}


#' Run RPSBLAST to generate domain architectures for proteins of interest
#'
#' @param rpsblast_path
#' @param db_search_path Path to the BLAST databases
#' @param db
#' @param query
#' @param evalue
#' @param out
#' @param num_threads
#'
#' @return
#' @export
#'
#' @examples
run_rpsblast <- function(rpsblast_path, db_search_path,
                         db = "refseq", query, evalue = "1e-5",
                         out, num_threads = 1) {
  #' Run RPSBLAST to generate domain architectures for proteins of interest
  #' @author Samuel Chen, Janani Ravi
  #' @param path
  #' @param db_search_path Path to the BLAST databases
  #' @param db
  #' @param query
  #' @param evalue
  #' @param out
  #' @param num_threads
  start <- Sys.time()
  system(paste0("export BLASTDB=/", db_search_path))
  system2(
    command = rpsblast_path,
    args = c(
      "-db", db,
      "-query", query,
      "-evalue", evalue,
      "-out", out,
      "-num_threads", num_threads
      #                  , "-outfmt", outfmt
    )
    print(Sys.time() - start)
}
