# Author(s): Samuel Chen
# Last modified: 2020

setClass("MolEvolData",
         slots = list(
             df = "data.frame",
             fasta_path = "character",
             msa_path = "character",
             ipr_path = "character",
             queries = "character",
             cln_path = "character",
             domainSeqs = "character",
             domain_ipr_path = "character"
         )
)
setClass("queryData",
         slots = list(
             df = "data.frame",
             fasta_path = "character",
             msa_path = "character",
             ipr_path = "character",
             queries = "character",
             cln_path = "character"
         )
)
setClass("blastUpload",
         slots = list(
             df = "character",
             seqs = "character"
         )
)
setClass("iprUpload",
         slots = list(
             df = "character",
             seqs = "character"
         )
)
setClass("seqUpload",
         slots = list(
             seqs = "character"
         )
)

combineFilesNopmap <- function(inpath, pattern, outpath,
                                 delim = "\t", skip = 0,
                                 col_names) {
    source_files <- dir(path = inpath, pattern = pattern, recursive = T)

    source_files_path <- paste0(inpath, source_files)

    dt_list <- map(source_files_path, function(x) {
        if (!grepl("query", x)) {
            fread(x, sep = delim, skip = skip, fill = T) # , col.names = col_names)
        }
    })

    combined <- rbindlist(dt_list, fill = T)

    fwrite(combined, outpath, sep = "\t")

    return(combined)
}


full_analysis_colnames <- c("AccNum", "QueryName", "STitle", "Species.x",
                            "TaxID.x", "Lineage", "PcPositive", "PcIdentity",
                            "AlnLength", "SAccNum", "SAllSeqID", "Mismatch",
                            "GapOpen", "QStart", "QEnd", "QLength", "SStart",
                            "SEnd", "SLength", "EValue",
                            "BitScore", "PcPosOrig",
                            "Name", "ClusterID", "DomArch.Pfam",
                            "DomArch.SMART",
                            "DomArch.CDD", "DomArch.TIGRFAM",
                            "DomArch.Phobius",
                            "DomArch.Gene3D", "DomArch.TMHMM",
                            "DomArch.SignalP_EUK",
                            "DomArch.SignalP_GRAM_NEGATIVE",
                            "DomArch.SignalP_GRAM_POSITIVE",
                            "Description", "Length", "TaxID.y",
                            "Species.y", "SourceDB", "Completeness")


processWrapperDir <- function(path, pinName, type = "full") {
    if (type == "full") {
        if (file.exists(paste0(path, "/cln_combined.tsv")) && file.exists(paste0(path, "/ipr_combined.tsv"))) {
            query_data <- read_tsv(paste0(path, "/query_data/query_data.full_analysis.tsv"))
            query_ipr_path <- paste0(path, "/query_data/query_data.iprscan_cln.tsv")
            query_msa_path <- paste0(path, "/query_data/query_seqs.msa")
            query_sequence_path <- paste0(path, "/query_data/query_data.all_accnums.fa")
            query_data$Query <- query_data$AccNum
            queries <- unique(query_data$Query)
            cln_path <- paste0(path, "/query_data/query_data.full_analysis.tsv")
            query_data <- query_data %>%
                mutate(QueryName = Name) %>%
                select(QueryName, Species, contains("Lineage"), contains("DomArch"), Query, Name) %>%
                distinct(QueryName, .keep_all = TRUE)
            query_wrapper_data <- new("queryData",
                                      df = query_data, fasta_path = query_sequence_path,
                                      msa_path = query_msa_path, ipr_path = query_ipr_path, queries = queries, cln_path = cln_path
            )
            com_blast_data <- read_tsv(path, "/blast_combined.tsv") %>% arrange(desc(PcPositive))
            ipr_blast_path <- paste0(path, "/ipr_combined.tsv")
            queries <- unique(com_blast_data$Query)
            com_blast_data$Lineage <- as.factor(com_blast_data$Lineage)
            com_blast_data <- com_blast_data %>% select(QueryName, everything())
            wrapper_data <- new("MolEvolData",
                                df = com_blast_data, queries = queries, ipr_path = ipr_blast_path, cln_path = com_blast_path,
                                fasta_path = query_sequence_path, msa_path = query_msa_path, domainSeqs = "", domain_ipr_path = ""
            )

            return(list(wrapper_data, query_wrapper_data))
        } else {

            # Get the query data
            query_data <- read_tsv(paste0(path, "/query_data/query_data.full_analysis.tsv"))
            query_ipr_path <- paste0(path, "/query_data/query_data.iprscan_cln.tsv")
            query_msa_path <- paste0(path, "/query_data/query_seqs.msa")
            query_sequence_path <- paste0(path, "/query_data/query_data.all_accnums.fa")
            query_data$Query <- query_data$AccNum
            # alignFasta(query_sequence_path, tool = "ClustalO", outpath = query_msa_path)
            queries <- unique(query_data$Query)
            cln_path <- paste0(path, "/query_data/query_data.full_analysis.tsv")
            query_data <- query_data %>%
                mutate(QueryName = Name) %>%
                select(QueryName, Species, contains("Lineage"), contains("DomArch"), Query, Name) %>%
                distinct(QueryName, .keep_all = TRUE)
            query_wrapper_data <- new("queryData",
                                      df = query_data, fasta_path = query_sequence_path,
                                      msa_path = query_msa_path, ipr_path = query_ipr_path, queries = queries, cln_path = cln_path
            )

            query_data <- query_data %>% select(Query, QueryName)

            com_blast_path <- paste0(path, "/blast_combined.tsv")
            ipr_blast_path <- paste0(path, "/ipr_combined.tsv")
            com_blast_data <- combineFilesNopmap(paste0(path, "/"),
                                                   pattern = "*.full_analysis.tsv", skip = 0,
                                                   col_names = c(), outpath = com_blast_path, delim = "\t"
            )

            ipr_blast_data <- combineFilesNopmap(paste0(path, "/"),
                                                   pattern = "*.iprscan_cln.tsv", skip = 0,
                                                   col_names = ipr_colnames, outpath = ipr_blast_path, delim = "\t"
            )
            com_blast_data <- merge(com_blast_data, query_data, by = "Query")
            com_blast_data <- com_blast_data %>% arrange(desc(PcPositive))
            fwrite(com_blast_data, com_blast_path, sep = "\t")

            queries <- unique(com_blast_data$Query)

            com_blast_data$Lineage <- as.factor(com_blast_data$Lineage)
            com_blast_data <- com_blast_data %>% select(QueryName, everything())

            wrapper_data <- new("MolEvolData",
                                df = com_blast_data, queries = queries, ipr_path = ipr_blast_path, cln_path = com_blast_path,
                                fasta_path = query_sequence_path, msa_path = query_msa_path, domainSeqs = "", domain_ipr_path = ""
            )

            return(list(wrapper_data, query_wrapper_data))
        }
    } else if (type == "dblast") {
        query_data <- read_tsv(paste0(path, "/query_data/query_data.full_analysis.tsv"))
        query_ipr_path <- paste0(path, "/query_data/query_data.iprscan_cln.tsv")
        query_msa_path <- paste0(path, "/query_data/query_seqs.msa")
        query_sequence_path <- paste0(path, "/query_data/query_data.all_accnums.fa")
        alignFasta(query_sequence_path, tool = "ClustalO", outpath = query_msa_path)
        cln_path <- paste0(path, "/query_data/query_data.full_analysis.tsv")
        query_data <- query_data %>% mutate(Query = AccNum)
        query_data <- query_data %>%
            mutate(QueryName = Name) %>%
            select(QueryName, Species, contains("Lineage"), contains("DomArch"),
                   Query, Name, AccNum) %>%
            distinct(QueryName, .keep_all = TRUE)
        queries <- unique(query_data$Query)
        query_wrapper_data <- new("queryData",
                                  df = query_data,
                                  fasta_path = query_sequence_path,
                                  msa_path = query_msa_path,
                                  ipr_path = query_ipr_path,
                                  queries = queries, cln_path = cln_path
        )
        query_data <- query_data %>% select(Query, QueryName)
        com_blast_path <- paste0(path, "/blast_combined.tsv")
        com_blast_data <- combineFilesNopmap(paste0(path, "/"),
                                               pattern = "*.blast.cln.tsv",
                                               skip = 0,
                                               col_names = c(),
                                               outpath = com_blast_path,
                                               delim = "\t"
        )
        com_blast_data <- merge(com_blast_data, query_data, by = "Query")
        queries <- unique(com_blast_data$Query)
        com_blast_data$Lineage <- as.factor(com_blast_data$Lineage)
        com_blast_data <- com_blast_data %>%
            select(QueryName, everything()) %>%
            arrange(desc(PcPositive))
        wrapper_data <- new("MolEvolData",
                            df = com_blast_data, queries = queries,
                            ipr_path = "", cln_path = com_blast_path,
                            fasta_path = query_sequence_path,
                            msa_path = query_msa_path, domainSeqs = "",
                            domain_ipr_path = ""
        )
        return(list(wrapper_data, query_wrapper_data))
    } else if (type == "phylo") {
        query_data <- read_tsv(paste0(path, "/query_data/query_data.full_analysis.tsv"))
        query_ipr_path <- paste0(path, "/query_data/query_data.iprscan_cln.tsv")
        query_msa_path <- paste0(path, "/query_data/query_seqs.msa")
        cln_path <- paste0(path, "/query_data/query_data.full_analysis.tsv")
        query_data <- query_data %>% mutate(Query = AccNum)
        query_data <- query_data %>%
            mutate(QueryName = Name) %>%
            select(QueryName, Species, contains("Lineage"), contains("DomArch"),
                   Query, Name, AccNum) %>%
            distinct(QueryName, .keep_all = TRUE)
        queries <- unique(query_data$Query)
        query_wrapper_data <- new("queryData",
                                  df = query_data,
                                  fasta_path = query_sequence_path,
                                  msa_path = query_msa_path,
                                  ipr_path = query_ipr_path,
                                  queries = queries, cln_path = cln_path
        )
        wrapper_data <- new("MolEvolData",
                            df = query_data,
                            fasta_path = query_sequence_path,
                            msa_path = query_msa_path,
                            ipr_path = query_ipr_path,
                            queries = queries, cln_path = cln_path,
                            domainSeqs = "", domain_ipr_path = ""
        )
        return(list(wrapper_data, query_wrapper_data))
    } else if (type == "da") {
        query_data <- read_tsv(paste0(path, "/query_data/query_data.full_analysis.tsv"))
        query_ipr_path <- paste0(path, "/query_data/query_data.iprscan_cln.tsv")
        query_msa_path <- paste0(path, "/query_data/query_seqs.msa")
        query_sequence_path <- paste0(path, "/query_data/query_data.all_accnums.fa")
        alignFasta(query_sequence_path, tool = "ClustalO", outpath = query_msa_path)
        cln_path <- paste0(path, "/query_data/query_data.full_analysis.tsv")
        query_data <- query_data %>%
            mutate(QueryName = Name) %>%
            select(QueryName, Species, contains("Lineage"), contains("DomArch"),
                   Query, Name) %>%
            distinct(QueryName, .keep_all = TRUE)
        queries <- unique(query_data$Query)
        query_wrapper_data <- new("queryData",
                                  df = query_data,
                                  fasta_path = query_sequence_path,
                                  msa_path = query_msa_path,
                                  ipr_path = query_ipr_path,
                                  queries = queries, cln_path = cln_path
        )
        com_blast_path <- paste0(path, "/blast_combined.tsv")
        com_blast_data <- combineFilesNopmap(paste0(path, "/"),
                                               pattern = "*.full_analysis.tsv",
                                               skip = 0,
                                               col_names = c(),
                                               outpath = com_blast_path,
                                               delim = "\t"
        )
        com_blast_data <- merge(com_blast_data, query_data, by = "Query")
        com_blast_data <- com_blast_data %>%
            select(QueryName, everything()) %>%
            arrange(desc(PcPositive))
        wrapper_data <- new("MolEvolData",
                            df = com_blast_data, queries = queries,
                            ipr_path = "", cln_path = com_blast_path,
                            fasta_path = query_sequence_path,
                            msa_path = query_msa_path,
                            domainSeqs = "", domain_ipr_path = ""
        )
        return(list(wrapper_data, query_wrapper_data))
    } else {
        stop("Unrecognized type. Please use 'full', 'dblast', 'phylo', or 'da'.")
    }
}

drop_empty <- function(df) {
    for (c in colnames(df)) {
        if (all(is.na(df[, ..c])) || all(df[, ..c] == "")) {
            df <- df %>% select(-contains(c))
        }
    }
    return(df)
}


clean_fetched <- function(df) {
    df <- as.data.frame(df)
    # df = df %>% cleanup_species()

    # Cleanup domarchs
    cols <- colnames(df)
    for (c in cols) {
        if (grepl("^DomArch", c)) {
            # Repeats
            old <- c
            new <- paste0(c, ".repeats")
            df <- df %>% cleanup_domarch(
                old = old, new = new,
                domains_rename = NULL, # domains_rename,
                domains_keep = NULL, # filter applied to only ClustName for now.
                domains_ignore = NULL, # !! should check and remove soon!
                repeat2s = FALSE,
                remove_tails = F, # new! check below if it works!
                remove_empty = F
            )

            new <- c
            df <- df %>% cleanup_domarch(
                old = old, new = new,
                domains_rename = NULL, # domains_rename,
                domains_keep = NULL, # filter applied to only ClustName for now.
                domains_ignore = NULL, # !! should check and remove soon!
                repeat2s = TRUE,
                remove_tails = F, # new! check below if it works!
                remove_empty = F
            )
        }
    }


    return(df)
}

# Description
# validation functions for accession numbers/headers for FASTA (validateAccNumFasta)
# & accession number input (validateAccNum) submission types.

library(Biostrings)
library(httr)
library(httr2)
library(rentrez)
library(shiny)

getAccNumFromFasta <- function(biostrings_aa_string_set, verbose = FALSE) {
    # parsing/cleaning accession numbers using the same methods as
    # `upstream_scripts/00_submit_full.R`'s `get_sequences()` function
    accnums <- c()
    for (header in names(biostrings_aa_string_set)) {
        # case 1 UnitProtKB formatted FASTA
        # header: https://www.uniprot.org/help/fasta-headers
        if (grepl("\\|", header)) {
            accnum <- unlist(strsplit(header, "\\|"))[2]
            # case 2 NCBI formatted FASTA
            # header: https://www.ncbi.nlm.nih.gov/genbank/fastaformat
        } else if (grepl(" ", header)) {
            accnum <- unlist(strsplit(header, " "))[1]
            # case 3 neither delimiter present; use the whole header
        } else {
            accnum <- header
        }
        # print debug info
        if (verbose == TRUE) {
            cat("header:", header, "\n")
            cat("accnum_parsed:", accnum, "\n")
        }
        accnums <- append(accnums, accnum)
    }
    print(accnums)
    return(accnums)
}

validateAccNumFasta <- function(text) {
    # INPUT: text from input object that houses FASTA data
    # validate the headers/accnums for FASTA submission
    # Return:
    #   T: when no duplicate accnums after parsing
    #   F: when duplicate accnums are present after parsing

    # TRY write tmp file for biostrings parsing
    path <- tryCatch(
        expr = {
            path <- tempfile()
            write(text, path)
            path
        },
        error = function(e) {
            NULL
        }
    )
    if (is.null(path)) {
        cat("failed to write temp fasta during accnum validation\n", file = stderr())
        return(FALSE) # validation fail
    }

    # TRY read seq
    fasta <- tryCatch(
        expr = {
            fasta <- Biostrings::readAAStringSet(path)
            fasta
        },
        error = function(e) {
            NULL
        },
        finally = {
            unlink(path)
        }
    )
    if (is.null(fasta)) {
        cat("failed to read sequence during accnum validation\n", file = stderr())
        return(FALSE) # validation fail
    }

    # parse accnums from fasta
    accnums <- getAccNumFromFasta(fasta)
    tb_accnum_counts <- tibble("frequencies" = table(accnums))
    # check for duplicates
    if (any(tb_accnum_counts$frequencies > 1)) {
        cat("duplicate headers found during fasta validation\n", file = stderr())
        return(FALSE) # validation fail
    } else {
        return(TRUE)
    }
}


# Test whether a single accession returns a valid protein from Entrez
isAccNumValidEntrez <- function(accnum, verbose = FALSE) {
    # empty accnum wil not raise an error from efetch, so test for this first
    if (nchar(accnum) <= 0) {if (verbose) {warning("empty accnum")}; return(FALSE)}

    # try performing a POST of the accnum, followed by a GET request for a protein
    # sequence
    # rentrez::entrez_fetch() returns an error when there's no protein,
    # so we're using try statements to handle this and return FALSE upon error
    result <- tryCatch(
        expr = {
            result <- rentrez::entrez_fetch(db = "protein", id = accnum, rettype = "fasta")
            if (verbose) {print(result)}
            TRUE
        },
        error = function(e) {if (verbose) {print(e)}; FALSE}
    )
    return(result)
}

# Test whether a single accession returns a valid protein from EBI
isAccNumValidEbi <- function(accnum, verbose = FALSE) {
    # validation: ensure there's some text to parse
    if (nchar(accnum) <= 0) {if (verbose) {warning("empty accnum")}; return(FALSE)}
    # construct a httr2 request to POST an accession number, then GET the fasta
    url_base <- "https://www.ebi.ac.uk/proteins/api/proteins?accession="
    url_protein <- paste0(url_base, accnum)
    req <- httr2::request(url_protein) |>
        httr2::req_headers("Accept" = "text/x-fasta", "Content-type" = "application/x-www-form-urlencoded")

    # wrap in try since a failed HTTP request will raise an R error
    try(httr2::req_perform(req))
    # get the HTTP response code
    resp <- httr2::last_response()
    # assign result based on HTTP response code ('200' is a success)
    result <- ifelse(resp$status == 200, TRUE, FALSE)
    if (verbose) {
        msg <- stringr::str_glue(
            "EBI protein query results:\n",
            "\taccnum: {accnum}\t validation result: {result}\n"
        ) |> print()
    }
    return(result)
}

# Perform a series of API reqs using NCBI entrez to validate accession numbers
performEntrezReqs <- function(accnums, verbose = FALSE, track_progress = FALSE) {
    # API guidelines docs
    # ebi: https://www.ebi.ac.uk/proteins/api/doc/index.html
    # entrez recommends no more than 3 POSTs per second
    # simple method, sleep for 1 second after every 3rd POST
    i <- 0
    results <- vapply(
        X = accnums,
        FUN = function(accnum) {
            result <- isAccNumValidEntrez(accnum, verbose = TRUE)
            i <<- i + 1L
            if (i >= 3L) {Sys.sleep(1); print('sleeping for entrez API reqs'); i <<- 0L}
            if (track_progress) {incProgress(1)}
            result
        },
        FUN.VALUE = logical(1)
    )
    return(results)
}

performEbiReqs <- function(accnums, verbose = FALSE, track_progress = FALSE) {
    # API guidelines docs
    # ebi: https://www.ebi.ac.uk/proteins/api/doc/index.html
    # EBI allows 200 POSTs per second
    # simple method, sleep for 1 second after every 200th POST
    i <- 0
    results <- vapply(
        X = accnums,
        FUN = function(accnum) {
            result <- isAccNumValidEbi(accnum, verbose = TRUE)
            i <<- i + 1L
            if (i >= 200L) {Sys.sleep(1); i <<- 0L}
            if (track_progress) incProgress(1)
            result
        },
        FUN.VALUE = logical(1)
    )
    return(results)
}

# Validate accession numbers from MolEvolvR user input
validateAccNum <- function(text, verbose = FALSE, track_progress = FALSE, n_steps = integer()) {
    # API guidelines docs
    # entrez https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
    # ebi: https://www.ebi.ac.uk/proteins/api/doc/index.html

    # get accnums from three possible delimiters
    # 1. comma separated values with any amount of space between them
    # 2. one or more space delimiter
    # 3. "\n" (newline) delimiter
    accnums <- unlist(strsplit(text, "\\s*,\\s*|\\s+|\\n|,"))

    warning()

    if (track_progress) {incProgress(1 / n_steps,
        message = "Testing for duplicate accession numbers . . .")}
    # fail fast for duplicates
    if (any(duplicated(accnums))) {
        warning("duplicate accesion numbers found")
        if (track_progress) {setProgress(n_steps,
            message = "Duplicate accession numbers found")}
        return(FALSE)
    }
    if (track_progress) {incProgress(1 / n_steps,
        message = "Please wait (searching NCBI's protein database) . . .")}
    # entrez API
    entrez_results <- performEntrezReqs(accnums, verbose = verbose)
    if (verbose) {
        stringr::str_glue("rentrez
                results:\n\t{paste0(entrez_results, collapse=",")}\n") |>
            print()
    }

    # EBI API
    # if all of entrez resulted in success, then skip EBI. Else, try EBI
    if (!all(entrez_results)) {
        ebi_results <- performEbiReqs(accnums, verbose = verbose)
        if (verbose) {
            stringr::str_glue("ebi results:\n\t{paste0(ebi_results,
                              collapse=",")}\n") |>
                print()
        }
        if (track_progress) {incProgress(1 / n_steps, message = "Please wait
            (searching EBI's protein database) . . .")}
        # OR test each accnum result across both dbs
        final_results <- ebi_results | entrez_results
    } else {
        if (track_progress) {setProgress(n_steps, message = "All accession
                                         numbers validated")}
        final_results <- entrez_results
    }
    failed_accnums <- accnums[which(!final_results)]
    if (length(failed_accnums) >= 1) {
        msg <- stringr::str_glue(
            "The following accnums failed validation:\n",
            "\t{paste0(failed_accnums, collapse=",")}\n"
        )
        warning(msg)
        if (track_progress) {setProgress(n_steps, message = "Accession number
                                         validation failed")}
    }
    setProgress(n_steps, message = "Accession number validation finished")
    return(
        list(
            "validation_result" = all(final_results),
            "failed_accnums" = failed_accnums
        )
    )
}

#===============================================================================
# Description
#   validation function for evalue input
#===============================================================================
validateEvalue <- function(input_value) {
    is_valid_evalue <- is.numeric(input_value) && input_value != 0
    return(is_valid_evalue)
}

#===============================================================================
# Author(s): JK
# Last modified: 2023_06
#===============================================================================
#-------------------------------------------------------------------------------
guessSeqType <- function(single_fasta, dna_guess_cutoff = 0.9, other_guess_cutoff = 0.5) {
    tb <- as_tibble(alphabetFrequency(single_fasta))
    n_other <- if ("other" %in% colnames(tb)) sum(unlist(tb["other"])) else 0

    aa_cols <- setdiff(Biostrings::AA_ALPHABET, c("*", "."))
    tb_dna <- tb %>% select(all_of(Biostrings::DNA_ALPHABET))
    tb_aa <- tb %>% select(all_of(aa_cols))

    total <- nchar(single_fasta)
    est_dna_prop <- sum(unlist(tb_dna)) / total
    est_aa_prop <- sum(unlist(tb_aa)) / total
    other_prop <- n_other / total
    cat(
        names(single_fasta), "\n",
        "estimated DNA alphabet proportion:", est_dna_prop, "\n",
        "estimated AA alphabet proportion:", est_aa_prop, "\n",
        "'other' alphabet proportion:", other_prop, "\n"
    )

    if (est_dna_prop >= dna_guess_cutoff) {
        guess <- "DNA"
    } else if (other_prop >= other_guess_cutoff) {
        guess <- NA
    } else {
        guess <- "AA"
    }
    cat("Guess:", guess, "\n")
    return(guess)
}

validateSeqBody <- function(text) {
    # convert string to single letter character vector
    individual_chars <- unlist(strsplit(text, ""))
    # return the characters from input that are NOT in the AA_ALPHABET
    invalid_chars <- setdiff(individual_chars, Biostrings::AA_ALPHABET)
    n_invalid_chars <- length(invalid_chars)

    # 0 length character vector means the sequence contains
    # characters that are all present in the AA_ALPHABET
    # >0 length character vector means there are some characters
    # which are NOT found in the AA_ALPHABET
    if (n_invalid_chars == 0) {
        return(TRUE)
    } else {
        warning("A sequence contains characters not found in the AA_ALPHABET")
        return(FALSE)
    }
}

#-------------------------------------------------------------------------------
validateFasta <- function(fasta_path, .type = "AA") {
    # handle case of fasta being unreadable/unrecognized format
    fasta <- tryCatch(
        expr = Biostrings::readAAStringSet(fasta_path),
        error = function(e) {
            warning(paste0("error: failed to read sequence from: ", fasta_path))
            return(NULL)
        }
    )
    # failure to read yields quick return
    if (is.null(fasta)) {
        return(FALSE) # validation fail
    }

    # require at least some characters for the sequence
    if (!(sum(nchar(unlist(fasta))) >= 1)) {
        warning(paste0("fasta does not have any sequence characters"))
        return(FALSE) # validation fail
    }
    print(fasta)

    ### validate sequence BODY
    # iteratively return boolean value for valid seq body (FALSE = invalid)
    seq_body_validations <- c()
    for (i in seq_along(1:length(fasta))) {
        seq_body <- as.character(fasta[i])
        seq_body_validations <- c(validateSeqBody(seq_body), seq_body_validations)
    }
    # require all sequences to have chars only in AA_ALPHABET
    if (!(all(seq_body_validations))) {
        warning("At least one sequence has non-IUPAC AA characters")
        return(FALSE) # validation fail
    }

    ### validate sequence TYPE
    # iteratively guess the type (DNA, AA, or other of each seq in the FASTA file)
    seq_types <- c()
    for (idx in seq_along(1:length(fasta))) {
        seq_types <- c(guessSeqType(fasta[idx]), seq_types)
    }

    switch(.type,
           "DNA" = {
               anti_type <- "AA"
           },
           "AA" = {
               anti_type <- "DNA"
           }
    )

    # if all seqs are of desired .type, only then TRUE
    if (sum(is.na(seq_types)) != 0) {
        return(FALSE) # validation fail
    } else if (anti_type %in% seq_types) {
        return(FALSE) # validation fail
    } else {
        return(TRUE) # validation success
    }
}

# Author(s): Samuel Chen, Lo Sosinski, Joe Burke
# Last modified: 2021

alpha_numeric <- c(0:9, letters, LETTERS)

rand_string <- function(length, characters = alpha_numeric, post_fix = "", ignorelist = c()) {
    str <- NULL

    while (is.null(str) || (str %in% ignorelist)) {
        chars <- sample(characters, size = length, replace = TRUE)
        str <- paste0(paste(chars, collapse = ""), post_fix)
    }

    return(str)
}
strsplit_vect <- function(strings, pattern = "", pos) {
    split_strings <- map(strsplit(strings, "_"), function(x) x[1]) %>% unlist()

    return(split_strings)
}
