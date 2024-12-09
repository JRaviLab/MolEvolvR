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
#' @importFrom rlang warn abort inform
#'
#' @return No return, but run to close all outstanding `sink()`s
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sinkReset()
#' }
sinkReset <- function() {
  # Handle all errors and warnings
  tryCatch({
    for (i in seq_len(sink.number())) {
      sink(NULL)
    }
    rlang::inform("All sinks closed", class = "sink_reset_info")
  }, error = function(e) {
    rlang::abort(paste("Error: ", e$message), class = "sink_reset_error")
  }, warning = function(w) {
    rlang::warn(paste("Warning: ", w$message), class = "sink_reset_warning")
  }, finally = {
    # If any additional cleanup is needed, it can be done here
    if (sink.number() > 0) {
      # Additional cleanup if sinks are still open
      rlang::inform("Some sinks remain open, ensure proper cleanup.",
                    class = "sink_cleanup_warning")
    }
  })
}


#' addLineage
#'
#' @param df A `data.frame` containing the input data. One column must contain 
#' the accession numbers.
#' @param acc_col A string specifying the column name in `df` that holds the 
#' accession numbers. Defaults to `"AccNum"`.
#' @param assembly_path A string specifying the path to the `assembly_summary.txt` 
#' file. This file contains metadata about assemblies.
#' @param lineagelookup_path A string specifying the path to the lineage lookup 
#' file, which contains a mapping from tax IDs to their corresponding lineages.
#' @param ipgout_path (Optional) A string specifying the path where IPG database 
#' fetch results will be saved. If `NULL`, the results are not written to a file.
#' @param plan A string specifying the parallelization strategy for the future
#' package, such as `"sequential"` or `"multisession"`.
#'
#' @importFrom dplyr pull
#' @importFrom magrittr %>%
#' @importFrom rlang sym warn abort inform
#'
#' @return A `data.frame` that combines the original `df` with the lineage 
#' information.
#' @export
#'
#' @examples
#' \dontrun{
#' addLineage()
#' }
addLineage <- function(df, acc_col = "AccNum", assembly_path,
                       lineagelookup_path, ipgout_path = NULL,
                       plan = "sequential", ...) {
  # check for validate inputs
  if (!is.data.frame(df)) {
    rlang::abort("Input 'df' must be a data frame.", class = "input_error")
  }

  if (!acc_col %in% colnames(df)) {
    rlang::abort(paste("Column", acc_col,
                       "not found in data frame."), class = "column_error")
  }

  # Ensure paths are character strings
  if (!is.character(assembly_path) || !is.character(lineagelookup_path)) {
    rlang::abort("Both 'assembly_path' and
                 'lineagelookup_path' must be character strings.",
                 class = "path_type_error")
  }

  # Ensure paths exist
  if (!file.exists(assembly_path)) {
    rlang::abort(paste("Assembly file not found at:",
                       assembly_path), class = "file_not_found_error")
  }

  if (!file.exists(lineagelookup_path)) {
    rlang::abort(paste("Lineage lookup file not found at:",
                       lineagelookup_path), class = "file_not_found_error")
  }
  tryCatch({
    # Attempt to add lineages
    acc_col <- sym(acc_col)
    accessions <- df %>% pull(acc_col)
    lins <- acc2Lineage(
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
    rlang::abort(paste("Error during lineage addition:", e$message),
                 class = "lineage_addition_error")
  }, warning = function(w) {
    rlang::warn(paste("Warning during lineage addition:", w$message),
                class = "lineage_addition_warning")
  })

}


#' acc2Lineage
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description This function combines 'efetchIPG()' and 'IPG2Lineage()' to map a set
#' of protein accessions to their assembly (GCA_ID), tax ID, and lineage.
#'
#' @param accessions Character vector of protein accessions
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the \link[MolEvolvR]{downloadAssemblySummary} function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' @param ipgout_path Path to write the results of the efetch run of the accessions
#' on the ipg database. If NULL, the file will not be written. Defaults to NULL
#' @param plan A string specifying the parallelization strategy for the future
#' package, such as `"sequential"` or `"multisession"`.
#'
#' @importFrom rlang warn abort inform
#'
#' @return A `data.table` that contains the lineage information, mapping protein 
#' accessions to their tax IDs and lineages.
#' @export
#'
#' @examples
#' \dontrun{
#' acc2Lineage()
#' }
acc2Lineage <- function(accessions, assembly_path,
                        lineagelookup_path, ipgout_path = NULL,
                        plan = "sequential", ...) {
  tmp_ipg <- F
  if (is.null(ipgout_path)) {
    tmp_ipg <- TRUE
    ipgout_path <- tempfile("ipg", fileext = ".txt")
  }

  lins <- NULL
  tryCatch({
    # Attempt to fetch IPG
    efetchIPG(accessions, out_path = ipgout_path, plan)

    # Attempt to process IPG to lineages
    lins <- IPG2Lineage(accessions, ipgout_path,
                        assembly_path, lineagelookup_path)
  }, error = function(e) {
    rlang::abort(
      message = paste("An error occurred during IPG fetching
                      or lineage processing:", e$message),
      class = "lineage_processing_error",
      # adding additional context
      accessions = accessions,
      assembly_path = assembly_path,
      lineagelookup_path = lineagelookup_path,
      ipgout_path = ipgout_path,
      plan = plan
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning during IPG fetching
                      or lineage processing:", w$message),
      class = "lineage_processing_warning",
      accessions = accessions,
      assembly_path = assembly_path,
      lineagelookup_path = lineagelookup_path,
      ipgout_path = ipgout_path,
      plan = plan
    )
  }, finally = {
    # Cleanup: delete temporary IPG file if it was created
    if (tmp_ipg && file.exists(ipgout_path)) {
      unlink(ipgout_path)
    }
  })

  return(lins)
}


#' efetchIPG
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Perform efetch on the ipg database and write the results to out_path
#' @param accnums Character vector containing the accession numbers to query on
#' the ipg database
#' @param out_path Path to write the efetch results to
#' @param plan A string specifying the parallelization strategy for the future
#' package, such as `"sequential"` or `"multisession"`.
#'
#' @importFrom furrr future_map
#' @importFrom future plan
#' @importFrom rentrez entrez_fetch
#' @importFrom rlang warn abort inform
#'
#' @return No return value. The function writes the fetched results to `out_path`.
#' @export
#'
#' @examples
#' \dontrun{
#' efetchIPG()
#' }
efetchIPG <- function(accnums, out_path, plan = "sequential", ...) {
  # Argument validation
  if (!is.character(accnums) || length(accnums) == 0) {
    rlang::abort("Error: 'accnums' must be a non-empty character vector.",
                 class = "validation_error")
  }

  if (!is.character(out_path) || nchar(out_path) == 0) {
    rlang::abort("Error: 'out_path' must be a non-empty string.",
                 class = "validation_error")
  }

  if (!is.function(plan)) {
    rlang::abort("Error: 'plan' must be a valid plan function.",
                 class = "validation_error")
  }
  if (length(accnums) > 0) {
    partition <- function(in_data, groups) {
      # \\TODO This function should be defined outside of efetchIPG().
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
      rlang::abort(
        message = paste("An error occurred: ", e$message),
        class = "fetch_error",
        accnums = accnums,
        out_path = out_path,
        plan = plan
      )
    }, warning = function(w) {
      rlang::warn(
        message = paste("Warning: ", w$message),
        class = "fetch_warning",
        accnums = accnums,
        out_path = out_path,
        plan = plan
      )
    }, finally = {
      # Ensure the sink is closed in case of errors
      if (sink.number() > 0) sink(NULL)
    })
  }
}



#' IPG2Lineage
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Takes the resulting file
#'              of an efetch run on the ipg database and
#'
#' @param accessions Character vector of protein accessions
#' @param ipg_file Filepath to the file containing results of an efetch run on the
#' ipg database. The protein accession in 'accessions' should be contained in this
#' file
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the \link[MolEvolvR]{downloadAssemblySummary} function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "createLineageLookup()" function
#'
#' @importFrom data.table fread
#' @importFrom rlang warn abort inform
#'
#' @return A `data.table` with the lineage information for the provided protein 
#' accessions.
#' @export
#'
#' @examples
#' \dontrun{
#' IPG2Lineage()
#' }
#'
IPG2Lineage <- function(accessions, ipg_file,
                        assembly_path, lineagelookup_path, ...) {
  # Argument validation for accessions
  if (!is.character(accessions) || length(accessions) == 0) {
    rlang::abort("Input 'accessions' must be a non-empty
          character vector.", class = "validation_error")
  }

  # check for validate inputs
  if (!is.character(ipg_file)) {
    rlang::abort("Input 'ipg_file' must be a
          character string.", class = "validation_error")
  }

  # Ensure paths are character strings
  if (!is.character(assembly_path) || !is.character(lineagelookup_path)) {
    rlang::abort("Both 'assembly_path' and 'lineagelookup_path'
                 must be character strings.", class = "validation_error")
  }

  # Ensure paths exist
  if (!file.exists(assembly_path)) {
    rlang::abort(paste("Assembly file not found at:", assembly_path),
                 class = "file_error")
  }

  if (!file.exists(lineagelookup_path)) {
    rlang::abort(paste("Lineage lookup file not found at:", lineagelookup_path),
                 class = "file_error")
  }

  # Process the IPG file
  try({
    # Attempt to read the IPG file
    ipg_dt <- fread(ipg_file, sep = "\t", fill = T)

    # Filter the IPG data table to only include the accessions
    ipg_dt <- ipg_dt[Protein %in% accessions]

    # Rename the 'Assembly' column to 'GCA_ID'
    ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")

    # Convert the IPG data table to a lineage data table
    lins <- GCA2Lineage(prot_data = ipg_dt, assembly_path, lineagelookup_path)

    # Filter out rows with missing lineage information
    lins <- lins[!is.na(Lineage)] %>% unique()

    return(lins)
  }, error = function(e) {
    rlang::abort(
      message = paste("An error occurred: ", e$message),
      class = "processing_error",
      accessions = accessions,
      ipg_file = ipg_file,
      assembly_path = assembly_path,
      lineagelookup_path = lineagelookup_path
    )
  }, warning = function(w) {
    rlang::warn(
      message = paste("Warning: ", w$message),
      class = "processing_warning",
      accessions = accessions,
      ipg_file = ipg_file,
      assembly_path = assembly_path,
      lineagelookup_path = lineagelookup_path
    )
  })

}





# efetchIPG <- function(accnums, outpath)
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
