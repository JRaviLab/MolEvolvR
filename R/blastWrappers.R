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
    )
    print(Sys.time() - start)
}
