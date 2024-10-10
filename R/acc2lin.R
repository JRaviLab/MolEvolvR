# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(rentrez))
# suppressPackageStartupMessages(library(future))
# suppressPackageStartupMessages(library(furrr))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(biomartr))

# https://stackoverflow.com/questions/18730491/sink-does-not-release-file
#' Sink Reset
#'
#' @return No return, but run to close all outstanding `sink()`s
#'         and handles any errors or warnings that occur during the process.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sink.reset()
#' }
sink.reset <- function() {
  # Handle all errors and warnings
  tryCatch({
    for (i in seq_len(sink.number())) {
      sink(NULL)
    }
    print("All sinks closed")
  }, error = function(e) {
    print(paste("Error: ", e$message))
  }, warning = function(w) {
    print(paste("Warning: ", w$message))
  }, finally = {
    print("resetSink function execution completed.")
  })
}


#' Add Lineages
#'
#' @param df
#' @param acc_col
#' @param assembly_path
#' @param lineagelookup_path
#' @param ipgout_path
#' @param plan
#'
#' @importFrom dplyr pull
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' add_lins()
#' }
add_lins <- function(df, acc_col = "AccNum", assembly_path,
                     lineagelookup_path, ipgout_path = NULL,
                     plan = "sequential") {
  # check for validate inputs
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data frame.")
  }

  if (!acc_col %in% colnames(df)) {
    stop(paste("Column", acc_col, "not found in data frame."))
  }

  # Ensure paths are character strings
  if (!is.character(assembly_path) || !is.character(lineagelookup_path)) {
    stop("Both 'assembly_path' and 
         'lineagelookup_path' must be character strings.")
  }

  # Ensure paths exist
  if (!file.exists(assembly_path)) {
    stop(paste("Assembly file not found at:", assembly_path))
  }

  if (!file.exists(lineagelookup_path)) {
    stop(paste("Lineage lookup file not found at:", lineagelookup_path))
  }
    tryCatch({
      # Attempt to add lineages
      acc_col <- sym(acc_col)
      accessions <- df %>% pull(acc_col)
      lins <- acc2lin(
        accessions, assembly_path, lineagelookup_path, ipgout_path, plan
      )

      # Drop a lot of the unimportant columns for now? 
      # will make merging much easier
      lins <- lins[, c(
        "Strand", "Start", "Stop", "Nucleotide Accession", "Source",
        "Id", "Strain"
      ) := NULL]
      lins <- unique(lins)

      # dup <- lins %>% group_by(Protein) %>% 
      # summarize(count = n()) %>% filter(count > 1) %>%
      # pull(Protein)

      merged <- merge(df, lins, by.x = acc_col, by.y = "Protein", all.x = TRUE)
      return(merged)
    }, error = function(e) {
      print(paste("Error: ", e$message))
    }, warning = function(w) {
      print(paste("Warning: ", w$message))
    }, finally = {
      print("addLineages function execution completed.")
    })

}


#' acc2lin
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description This function combines 'efetch_ipg()'
#'              and 'ipg2lin()' to map a set
#' of protein accessions to their assembly (GCA_ID), tax ID, and lineage.
#'
#' @param accessions Character vector of protein accessions
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "DownloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' @param ipgout_path Path to write the results 
#'                    of the efetch run of the accessions
#' on the ipg database. If NULL, the file will not be written. Defaults to NULL
#' @param plan
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' acc2lin()
#' }
acc2lin <- function(accessions, assembly_path, 
                    lineagelookup_path, ipgout_path = NULL, 
                    plan = "sequential") {
  tmp_ipg <- F
  if (is.null(ipgout_path)) {
    tmp_ipg <- T
    ipgout_path <- tempfile("ipg", fileext = ".txt")
  }

  lins <- NULL
  tryCatch({
    # Attempt to fetch IPG
    efetch_ipg(accessions, out_path = ipgout_path, plan)

    # Attempt to process IPG to lineages
    lins <- ipg2lin(accessions, ipgout_path, assembly_path, lineagelookup_path)
  }, error = function(e) {
    print(paste("An error occurred: ", e$message))
  }, warning = function(w) {
    print(paste("Warning: ", w$message))
  }, finally = {
    print("acc2lin function execution completed.")
  })

  if (tmp_ipg) {
    unlink(tempdir(), recursive = T)
  }
  return(lins)
}


#' efetch_ipg
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Perform efetch on the ipg database
#'              and write the results to out_path
#'
#' @param accnums Character vector containing the accession numbers to query on
#' the ipg database
#' @param out_path Path to write the efetch results to
#' @param plan
#'
#' @importFrom furrr future_map
#' @importFrom future plan
#' @importFrom rentrez entrez_fetch
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' efetch_ipg()
#' }
efetch_ipg <- function(accnums, out_path, plan = "sequential") {
  # Argument validation
  if (!is.character(accnums) || length(accnums) == 0) {
    stop("Error: 'accnums' must be a non-empty character vector.")
  }

  if (!is.character(out_path) || nchar(out_path) == 0) {
    stop("Error: 'out_path' must be a non-empty string.")
  }

  if (!is.function(plan)) {
    stop("Error: 'plan' must be a valid plan function.")
  }
  if (length(accnums) > 0) {
    partition <- function(in_data, groups) {
      # \\TODO This function should be defined outside of efetch_ipg().
      # It can be non-exported/internal
      # Partition data to limit number of queries per second for rentrez fetch:
      # limit of 10/second w/ key
      l <- length(in_data)

      partitioned <- list()
      for (i in 1:groups){
        partitioned[[i]] <- in_data[seq.int(i, l, groups)]
      }

      return(partitioned)
    }
    tryCatch({
      # Set the future plan strategy
      plan(strategy = plan, .skip = T)


      min_groups <- length(accnums) / 200
      groups <- min(max(min_groups, 15), length(accnums))
      partitioned_acc <- partition(accnums, groups)

      # Open the sink to the output path
      sink(out_path)

      a <- future_map(1:length(partitioned_acc), function(x) {
        # Avoid hitting the rate API limit
        if (x %% 9 == 0) {
          Sys.sleep(1)
        }
        cat(
          entrez_fetch(
            id = partitioned_acc[[x]],
            db = "ipg",
            rettype = "xml",
            api_key = "YOUR_KEY_HERE" ## Can this be included in public package?
          )
        )
      })
      sink(NULL)
    }, error = function(e) {
      print(paste("An error occurred: ", e$message))
    }, warning = function(w) {
      print(paste("Warning: ", w$message))
    }, finally = {
      print("efetch_ipg function execution completed.")
    })
  }
}



#' ipg2lin
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Takes the resulting file
#'              of an efetch run on the ipg database and
#'
#' @param accessions Character vector of protein accessions
#' @param ipg_file Filepath to the file
#'                 containing results of an efetch run on the
#' ipg database. The protein accession in
#'               'accessions' should be contained in this
#' file
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "DownloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "createLineageLookup()" function
#'
#' @importFrom data.table fread
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' ipg2lin()
#' }
#'
ipg2lin <- function(accessions, ipg_file, assembly_path, lineagelookup_path) {
  # Argument validation for accessions
  if (!is.character(accessions) || length(accessions) == 0) {
    stop("Input 'accessions' must be a non-empty character vector.")
  }

  # check for validate inputs
  if (!is.character(ipg_file)) {
    stop("Input 'ipg_file' must be a character string.")
  }
  # Ensure paths are character strings
  if (!is.character(assembly_path) || !is.character(lineagelookup_path)) {
    stop("Both 'assembly_path' and 
         'lineagelookup_path' must be character strings.")
  }

  # Ensure paths exist
  if (!file.exists(assembly_path)) {
    stop(paste("Assembly file not found at:", assembly_path))
  }

  if (!file.exists(lineagelookup_path)) {
    stop(paste("Lineage lookup file not found at:", lineagelookup_path))
  }

  try({
    # Attempt to read the IPG file
    ipg_dt <- fread(ipg_file, sep = "\t", fill = T)

    # Filter the IPG data table to only include the accessions
    ipg_dt <- ipg_dt[Protein %in% accessions]

    # Rename the 'Assembly' column to 'GCA_ID'
    ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")

    # Convert the IPG data table to a lineage data table
    lins <- GCA2Lins(prot_data = ipg_dt, assembly_path, lineagelookup_path)

    # Filter out rows with missing lineage information
    lins <- lins[!is.na(Lineage)] %>% unique()

    return(lins)
  }, error = function(e) {
    print(paste("An error occurred: ", e$message))
  }, warning = function(w) {
    print(paste("Warning: ", w$message))
  }, finally = {
    print("ipg2lin function execution completed.")
  })
}





# efetch_ipg <- function(accnums, outpath)
# {
#   SIZE = 250
#   lower_bound = 1
#   groups = ceiling(length(accnums)/250)
#   for(i in 1:groups)
#   {
#     upper_bound = min((lower_bound+250),length(accnums))
#
#     sub_acc = accnums[lower_bound:upper_bound]
#
#     lower_bound = upper_bound +1
#
#     # Post needs things UID's, and doesn't really take AccNums
#     if(i == 1)
#     {
#       acc_webhistory <- entrez_post(db = "ipg", id =  sub_acc,api_key = "YOUR_KEY_HERE")
#     }
#     else
#     {
#       acc_webhistory <- entrez_post(db = "ipg", id =  sub_acc, web_history = acc_webhistory,
#                                     api_key = "YOUR_KEY_HERE")
#     }
#   }
#   browser()
#   sink(outpath)
#   cat(entrez_fetch("ipg", rettype = "xml", web_history =  acc_webhistory))
#   sink(NULL)
# }
