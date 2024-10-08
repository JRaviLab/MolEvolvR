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
#' @export
#'
#' @examples
#' \dontrun{
#' sinkReset()
#' }
sinkReset <- function() {
    for (i in seq_len(sink.number())) {
        sink(NULL)
    }
}


#' addLineage
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
#' addLineage()
#' }
addLineage <- function(df, acc_col = "AccNum", assembly_path,
    lineagelookup_path, ipgout_path = NULL, plan = "sequential", ...) {
    s_acc_col <- sym(acc_col)
    accessions <- df %>% pull(acc_col)
    lins <- acc2Lineage(accessions, assembly_path, lineagelookup_path, ipgout_path, plan)

    # Drop a lot of the unimportant columns for now? will make merging much easier
    lins <- lins[, c(
        "Strand", "Start", "Stop", "Nucleotide Accession", "Source",
        "Id", "Strain"
    ) := NULL]
    lins <- unique(lins)

    # dup <- lins %>% group_by(Protein) %>% summarize(count = n()) %>% filter(count > 1) %>%
    #   pull(Protein)

    merged <- merge(df, lins, by.x = acc_col, by.y = "Protein", all.x = TRUE)
    return(merged)
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
#' @param plan
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' acc2Lineage()
#' }
acc2Lineage <- function(accessions, assembly_path, lineagelookup_path, ipgout_path = NULL, plan = "sequential", ...) {
    tmp_ipg <- F
    if (is.null(ipgout_path)) {
        tmp_ipg <- T
        ipgout_path <- tempfile("ipg", fileext = ".txt")
    }
    efetchIPG(accessions, out_path = ipgout_path, plan)

    lins <- IPG2Lineage(accessions, ipgout_path, assembly_path, lineagelookup_path)

    if (tmp_ipg) {
        unlink(tempdir(), recursive = T)
    }
    return(lins)
}

#' efetchIPG
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Perform efetch on the ipg database and write the results to out_path
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
#' efetchIPG()
#' }
efetchIPG <- function(accnums, out_path, plan = "sequential", ...) {
    if (length(accnums) > 0) {
        partition <- function(in_data, groups) {
            # \\TODO This function should be defined outside of efetchIPG(). It can be non-exported/internal
            # Partition data to limit number of queries per second for rentrez fetch:
            # limit of 10/second w/ key
            l <- length(in_data)

            partitioned <- list()
            for (i in 1:groups)
            {
                partitioned[[i]] <- in_data[seq.int(i, l, groups)]
            }

            return(partitioned)
        }

        plan(strategy = plan, .skip = T)


        min_groups <- length(accnums) / 200
        groups <- min(max(min_groups, 15), length(accnums))
        partitioned_acc <- partition(accnums, groups)
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
    }
}

#' IPG2Lineage
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @description Takes the resulting file of an efetch run on the ipg database and
#'
#' @param accessions Character vector of protein accessions
#' @param ipg_file Filepath to the file containing results of an efetch run on the
#' ipg database. The protein accession in 'accessions' should be contained in this
#' file
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the \link[MolEvolvR]{downloadAssemblySummary} function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "create_lineage_lookup()" function
#'
#' @importFrom data.table fread
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' IPG2Lineage()
#' }
#'
IPG2Lineage <- function(accessions, ipg_file, assembly_path, lineagelookup_path, ...) {
    ipg_dt <- fread(ipg_file, sep = "\t", fill = T)

    ipg_dt <- ipg_dt[Protein %in% accessions]

    ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")

    lins <- GCA2Lins(prot_data = ipg_dt, assembly_path, lineagelookup_path)
    lins <- lins[!is.na(Lineage)] %>% unique()

    return(lins)
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
