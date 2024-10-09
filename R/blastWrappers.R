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

  # Argument validation
  if (!file.exists(deltablast_path)) {
    stop("The DELTABLAST executable path is invalid: ", deltablast_path)
  }
  if (!dir.exists(db_search_path)) {
    stop("The database search path is invalid: ", db_search_path)
  }
  if (!file.exists(query)) {
    stop("The query file path is invalid: ", query)
  }
  if (!is.numeric(as.numeric(evalue)) || as.numeric(evalue) <= 0) {
    stop("The evalue must be a positive number: ", evalue)
  }
  if (!is.numeric(num_alignments) || num_alignments <= 0) {
    stop("The number of alignments must be a 
         positive integer: ", num_alignments)
  }
  if (!is.numeric(num_threads) || num_threads <= 0) {
    stop("The number of threads must be a positive integer: ", num_threads)
  }

  start <- Sys.time()

  tryCatch({
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
  }, error = function(e) {
    message(paste("Error in run_deltablast: ", e))
  }, warning = function(w) {
    message(paste("Warning in run_deltablast: ", w))
  }, finally = {
    message("run_deltablast completed")
  })

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
  # Argument validation
  if (!file.exists(rpsblast_path)) {
    stop("The RPSBLAST executable path is invalid: ", rpsblast_path)
  }
  if (!dir.exists(db_search_path)) {
    stop("The database search path is invalid: ", db_search_path)
  }
  if (!file.exists(query)) {
    stop("The query file path is invalid: ", query)
  }
  if (!is.numeric(as.numeric(evalue)) || as.numeric(evalue) <= 0) {
    stop("The evalue must be a positive number: ", evalue)
  }
  if (!is.numeric(num_threads) || num_threads <= 0) {
    stop("The number of threads must be a positive integer: ", num_threads)
  }

  start <- Sys.time()

  tryCatch({

    system(paste0("export BLASTDB=/", db_search_path))

    system2(
      command = rpsblast_path,
      args = c(
        "-db", db,
        "-query", query,
        "-evalue", evalue,
        "-out", out,
        "-num_threads", num_threads
      )
    )
    print(Sys.time() - start)
  }, error = function(e) {
    message(paste("Error in run_rpsblast: ", e))
  }, warning = function(w) {
    message(paste("Warning in run_rpsblast: ", w))
  }, finally = {
    message("run_rpsblast completed")
  })

}
