# example data to test these functions is located at
# molevol_scripts/tests
# notes:
# - a protein with no domains (unlikely) found from
# interproscan CLI will return a completely empty file (0Bytes)

#' runIPRScan
#'
#' Run InterProScan on a given FASTA file and save the results to an
#' output file.
#'
#' @param filepath_fasta A string representing the path to the input FASTA file.
#' @param filepath_out A string representing the base path for the output file.
#' @param appl A character vector specifying the InterProScan applications to
#' use (e.g., "Pfam", "Gene3D"). Default is `c("Pfam", "Gene3D")`.
#'
#' @importFrom stringr str_glue
#'
#' @return A data frame containing the results from the InterProScan output
#' TSV file.
#'
#' @examples
#' \dontrun{
#' results <- runIPRScan(
#'     filepath_fasta = "path/to/your_fasta_file.fasta",
#'     filepath_out = "path/to/output_file",
#'     appl = c("Pfam", "Gene3D")
#' )
#' results
#' }
runIPRScan <- function(
        filepath_fasta,
        filepath_out, # do not inlucde file extension since ipr handles this
        appl = c("Pfam", "Gene3D")
        # destPartition = "LocalQ",
        # destQoS = "shortjobs"
    ) {
    # construct interproscan command
    cmd_iprscan <- stringr::str_glue(
        "iprscan -i {filepath_fasta} -b {filepath_out} --cpu 4 -f TSV ",
        "--appl {appl}"
    )
    # execute
    exit_code <- system(cmd_iprscan)
    if (exit_code != 0L) {
        warning("interproscan exited with non-zero code")
        return(NULL)
    }
    # read and return results
    df_iprscan <- readIPRScanTSV(paste0(filepath_out, ".tsv"))
    return(df_iprscan)
}

#' Constructor function for interproscan column names
#' (based upon the global variable written in
#' molevol_scripts/R/colnames_molevol.R)
#'
#' @return [chr] interproscan column names used throughout molevolvr
getIPRScanColNames <- function() {
    column_names <- c(
        "AccNum", "SeqMD5Digest", "SLength", "Analysis",
        "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
        "Status", "RunDate", "IPRAcc", "IPRDesc"
    )
    return(column_names)
}

#' construct column types for reading interproscan output TSVs
#' (based upon the global variable written in
#' molevol_scripts/R/colnames_molevol.R)
#' @return [collector] a named vector of type expecatations
#' for interproscan columns
#'
getIPRScanColTypes <- function() {
    column_types <- readr::cols(
        "AccNum" = readr::col_character(),
        "SeqMD5Digest" = readr::col_character(),
        "SLength" = readr::col_integer(),
        "Analysis" = readr::col_character(),
        "DB.ID" = readr::col_character(),
        "SignDesc" = readr::col_character(),
        "StartLoc" = readr::col_integer(),
        "StopLoc" = readr::col_integer(),
        "Score" = readr::col_double(),
        "Status" = readr::col_character(),
        "RunDate" = readr::col_character(),
        "IPRAcc" = readr::col_character(),
        "IPRDesc" = readr::col_character(),
    )
    return(column_types)
}

#' Read an interproscan output TSV with standardized
#' column names and types
#'
#' @param filepath [chr] path to interproscan output TSV
#'
#' @importFrom readr read_tsv
#'
#' @return [tbl_df] interproscan output table
readIPRScanTSV <- function(filepath) {
    df_ipr <- readr::read_tsv(filepath,
        col_types = getIPRScanColTypes(),
        col_names = getIPRScanColNames()
    )
    return(df_ipr)
}

#' For a given accession number, get the domain sequences using a interproscan
#' output table & the original FASTA file
#'
#' @param accnum [chr] a *single* accession number from the original fasta (fasta param)
#' which will be used to search for its sequence's domains (df_iprscan param)
#' @param fasta [AAStringSet] original fasta file which was fed into interproscan
#' @param df_iprscan [tbl_df] the output TSV of interproscan, read as a tibble with
#' readIPRScanTSV()
#' @param analysis [chr] the domain databases to extract sequences from
#'
#' @importFrom dplyr arrange filter mutate rowwise relocate select ungroup
#' @importFrom stringr str_glue
#' @importFrom XVector subseq
#'
#' @return [tbl_df] table with each domain sequence and a new identifier column
#'
#' @examples
#' \dontrun{
#' path_molevol_scripts <- file.path(Sys.getenv("DEV", unset = "/data/molevolvr_transfer/molevolvr_dev"), "molevol_scripts")
#' setwd(path_molevol_scripts)
#' source("R/fa2domain.R")
#' fasta <- Biostrings::readAAStringSet("./tests/example_protein.fa")
#' df_iprscan <- readIPRScanTSV("./tests/example_iprscan_valid.tsv")
#' accnum <- df_iprscan$AccNum[1]
#' df_iprscan_domains <- createIPRScanDomainTable(accnum, fasta, df_iprscan)
#' }
#'
createIPRScanDomainTable <- function(
        accnum,
        fasta,
        df_iprscan,
        analysis = c("Pfam", "Gene3D")) {
    cat("accnum:", accnum, "\n")
    # handle rare case of an interproscan result with no domains at all;
    # return the tibble with 0 rows quickly
    if (nrow(df_iprscan) < 1) {
        return(df_iprscan)
    }
    # filter for "Analysis" argument (i.e., the interproscan member database)
    # by default we're only selecting domains from a subset of the
    # interproscan member databases, but this param can be overwritten
    # and
    # filter for the accnum of interest (note: it's possible the accession
    # number is not in the table [i.e., it had no domains])
    df_iprscan_accnum <- df_iprscan |>
        dplyr::filter(.data$Analysis %in% analysis) |>
        dplyr::filter(.data$AccNum == accnum) |>
        dplyr::select(dplyr::all_of(c("AccNum", "DB.ID", "StartLoc", "StopLoc"))) |>
        dplyr::arrange(.data$StartLoc)
    # handle the case of no records after filtering by "Analysis"; return the tibble
    # with 0 rows quickly
    if (nrow(df_iprscan_accnum) < 1) {
        return(df_iprscan_accnum)
    }

    # create a new column to store the domain sequences
    df_iprscan_domains <- df_iprscan_accnum |>
        dplyr::rowwise() |>
        dplyr::mutate(
            seq_domain = XVector::subseq(
                fasta[[grep(pattern = .data$AccNum, x = names(fasta), fixed = TRUE)]],
                start = .data$StartLoc,
                end = .data$StopLoc
            ) |>
                as.character()
        )

    # create identifiers for each domain sequence
    df_iprscan_domains <- df_iprscan_domains |>
        dplyr::mutate(
            id_domain = stringr::str_glue("{AccNum}-{DB.ID}-{StartLoc}_{StopLoc}")
        ) |>
        dplyr::ungroup() |>
        dplyr::relocate(.data$id_domain, .before = 1)
    return(df_iprscan_domains)
}

#' Using the table returned from createIPRScanDomainTable, construct a
#' domain fasta for a single accession number in the original fasta
#' (i.e., the original fasta argument to createIPRScanDomainTable())
#'
#' @param df_iprscan_domains [tbl_df] return value from createIPRScanDomainTable
#'
#' @importFrom Biostrings AAStringSet
#' @importFrom dplyr mutate rowwise
#'
#' @return [AAStringSet] A domain fasta containing all the domains for a
#' single protein in the original fasta passed as an argument to createIPRScanDomainTable()
#'
#' @examples
#' \dontrun{
#' path_molevol_scripts <- file.path(Sys.getenv("DEV", unset = "/data/molevolvr_transfer/molevolvr_dev"), "molevol_scripts")
#' setwd(path_molevol_scripts)
#' source("R/fa2domain.R")
#' fasta <- Biostrings::readAAStringSet("./tests/example_protein.fa")
#' df_iprscan <- readIPRScanTSV("./tests/example_iprscan_valid.tsv")
#' accnum <- df_iprscan$AccNum[1]
#' df_iprscan_domains <- createIPRScanDomainTable(accnum, fasta, df_iprscan)
#' fasta_domains <- df_iprscan_domains |> convertIPRScanDomainTable2FA()
#' }
#'
convertIPRScanDomainTable2FA <- function(df_iprscan_domains) {
    # if there are no records (e.g., after filtering for Pfam analysis only)
    # the quickly return an empty AAStringSet object
    if (nrow(df_iprscan_domains) < 1) {
        return(Biostrings::AAStringSet())
    }

    # function with side effect to append a fasta when applied to tibble rows
    # 1. constructs a new fasta for a single domain and
    # 2. appends it to a `fasta_domains` object in the parent env
    # 3. returns the index of the new record; although,
    # note: return value is not particularly important here since we're
    # writing the main object within this function block
    fasta_domains <- Biostrings::AAStringSet()
    append_fasta_domains <- function(new_seq, new_seq_id) {
        idx_new_record <- length(fasta_domains) + 1
        fasta_domain <- Biostrings::AAStringSet(new_seq)
        names(fasta_domain) <- new_seq_id
        fasta_domains <<- c(fasta_domains, fasta_domain)
        return(idx_new_record)
    }

    # apply append_fasta_domains() rowwise
    df_iprscan_domains <- df_iprscan_domains |>
        dplyr::rowwise() |>
        dplyr::mutate(
            idx_new_record = append_fasta_domains(
                new_seq = seq_domain,
                new_seq_id = id_domain
            )
        )
    return(fasta_domains)
}

#' getDomainsFromFA
#'
#' @param fasta [AAStringSet] a protein (AA) fasta
#' @param df_iprscan [tbl_df] the interproscan results from the original fasta
#' @param analysis [chr] the domain databases to extract sequences from
#'
#' @importFrom Biostrings AAStringSet
#' @importFrom stringr str_glue
#'
#' @return fasta_domains [AAStringSet] fasta of domains
#'
#' @examples
#' \dontrun{
#' path_molevol_scripts <- file.path(Sys.getenv("DEV", unset = "/data/molevolvr_transfer/molevolvr_dev"), "molevol_scripts")
#' setwd(path_molevol_scripts)
#' source("R/fa2domain.R")
#' fasta <- Biostrings::readAAStringSet("./tests/example_protein.fa")
#' df_iprscan <- readIPRScanTSV("./tests/example_iprscan_valid.tsv")
#' getDomainsFromFA(fasta, df_iprscan)
#' }
#'
getDomainsFromFA <- function(
        fasta,
        df_iprscan,
        analysis = c("Pfam", "Gene3D"),
        verbose = FALSE) {
    # initialize an AAStringSet which will store
    # all the domain sequences
    # named "parent" since this will be appended
    # from a child environment
    parent_fasta_domains <- Biostrings::AAStringSet()

    # for each seq in the fasta use the interproscan
    # table construct a new fasta of the domains
    #
    # a boolean result keeps track of which sequences do or
    # do not have domain data that was added to the final
    # domain fasta
    results <- vapply(
        X = names(fasta),
        FUN = function(header) {
            # parse the accession number from header
            df_iprscan_domains <- createIPRScanDomainTable(
                header,
                fasta,
                df_iprscan,
                analysis = c("Pfam", "Gene3D")
            )
            # if the interpro results are empty OR
            # there's no domains for the analyses (databases)
            # then return early and do not append
            if (nrow(df_iprscan_domains) < 1) {
                if (verbose) {
                    msg <- stringr::str_glue(
                        "accession number: {header} had no domains for the ",
                        "selected analyes: {paste(analysis, collapse = ',')}\n"
                    )
                    warning(msg)
                }
                return(FALSE)
            }
            fasta_domains <- convertIPRScanDomainTable2FA(df_iprscan_domains)
            parent_fasta_domains <<- c(parent_fasta_domains, fasta_domains)
            return(TRUE)
        },
        FUN.VALUE = logical(1)
    )
    if (verbose) {
        msg <- stringr::str_glue(
            "{sum(results)} / {length(fasta)} accession numbers had ",
            "at least 1 domain available from the selected analyses: ",
            "{paste(analysis, collapse = ',')}\n"
        )
        print(msg)
    }
    return(parent_fasta_domains)
}
