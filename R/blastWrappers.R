# suppressPackageStartupMessages(library(tidyverse))

#' Run DELTABLAST to find homologs for proteins of interest
#'
#' @author Samuel Chen, Janani Ravi
#' @description
#' This function executes a Delta-BLAST search using the specified parameters
#' and database. It sets the BLAST database path, runs the Delta-BLAST command
#' with the given query, and outputs the results.
#'
#' @param deltablast_path Path to the Delta-BLAST executable.
#' @param db_search_path Path to the BLAST databases.
#' @param db Name of the BLAST database to search against (default is "refseq").
#' @param query Path to the input query file.
#' @param evalue E-value threshold for reporting matches (default is "1e-5").
#' @param out Path to the output file where results will be saved.
#' @param num_alignments Number of alignments to report.
#' @param num_threads Number of threads to use for the search (default is 1).
#'
#' @importFrom rlang warn abort inform
#'
#' @return This function does not return a value; it outputs results to the
#' specified file.
#' @export
#'
#' @examples
#' \dontrun{
#' runDeltaBlast(runDeltaBlast, db_search_path)
#' }
runDeltaBlast <- function(runDeltaBlast, db_search_path,
    db = "refseq", query, evalue = "1e-5",
    out, num_alignments, num_threads = 1) {

  # Argument validation
  if (!file.exists(deltablast_path)) {
    rlang::abort(paste("The DELTABLAST executable path is invalid:",
                       deltablast_path))
  }
  if (!dir.exists(db_search_path)) {
    rlang::abort(paste("The database search path is invalid:", db_search_path))
  }
  if (!file.exists(query)) {
    rlang::abort(paste("The query file path is invalid:", query))
  }
  if (!is.numeric(as.numeric(evalue)) || as.numeric(evalue) <= 0) {
    rlang::abort(paste("The evalue must be a positive number:", evalue))
  }
  if (!is.numeric(num_alignments) || num_alignments <= 0) {
    rlang::abort(paste("The number of alignments must be a positive integer:",
                       num_alignments))
  }
  if (!is.numeric(num_threads) || num_threads <= 0) {
    rlang::abort(paste("The number of threads must be a positive integer:",
                       num_threads))
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
    rlang::abort(
      message = paste("Error in run_deltablast:", e$message),
      class = "processing_error",
      deltablast_path = deltablast_path,
      db_search_path = db_search_path,
      query = query,
      out = out,
      num_alignments = num_alignments,
      num_threads = num_threads
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning in run_deltablast:", w$message),
      class = "processing_warning",
      deltablast_path = deltablast_path,
      db_search_path = db_search_path,
      query = query,
      out = out,
      num_alignments = num_alignments,
      num_threads = num_threads
    )
  })
}


#' Run RPSBLAST to generate domain architectures for proteins of interest
#'
#' @description
#' This function executes an RPS-BLAST search to generate domain architectures
#' for specified proteins. It sets the BLAST database path, runs the RPS-BLAST
#' command with the provided query, and outputs the results.
#'
#' @param rpsblast_path Path to the RPS-BLAST executable.
#' @param db_search_path Path to the BLAST databases.
#' @param db Name of the BLAST database to search against (default is "refseq").
#' @param query Path to the input query file.
#' @param evalue E-value threshold for reporting matches (default is "1e-5").
#' @param out Path to the output file where results will be saved.
#' @param num_threads Number of threads to use for the search (default is 1).
#'
#' @importFrom rlang warn abort inform
#'
#' @return This function does not return a value; it outputs results to the
#' specified file.
#' @export
#'
#' @examples
#' \dontrun{
#' runRSPBlast(rpsblast_path, db_search_path, query, out)
#' }
runRPSBlast <- function(rpsblast_path, db_search_path,
    db = "refseq", query, evalue = "1e-5",
    out, num_threads = 1) {
  
  # Argument validation
  if (!file.exists(rpsblast_path)) {
    rlang::abort(paste("The RPSBLAST executable path is invalid:",
                       rpsblast_path),
                 class = "file_error")
  }
  if (!dir.exists(db_search_path)) {
    rlang::abort(paste("The database search path is invalid:", db_search_path),
                 class = "file_error")
  }
  if (!file.exists(query)) {
    rlang::abort(paste("The query file path is invalid:", query),
                 class = "file_error")
  }
  if (!is.numeric(as.numeric(evalue)) || as.numeric(evalue) <= 0) {
    rlang::abort(paste("The evalue must be a positive number:", evalue),
          class = "validation_error")
  }
  if (!is.numeric(num_threads) || num_threads <= 0) {
    rlang::abort(paste("The number of threads must be a positive integer:",
                       num_threads),
                 class = "validation_error")
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
    rlang::abort(
      message = paste("Error in run_rpsblast:", e$message),
      class = "processing_error",
      rpsblast_path = rpsblast_path,
      db_search_path = db_search_path,
      query = query,
      out = out,
      num_threads = num_threads
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning in run_rpsblast:", w$message),
      class = "processing_warning",
      rpsblast_path = rpsblast_path,
      db_search_path = db_search_path,
      query = query,
      out = out,
      num_threads = num_threads
    )
  })

}
