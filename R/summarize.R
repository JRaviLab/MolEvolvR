## To summarize by lineages, DA and GC
## Created: Jun 07, 2019
## Modified: Dec 12, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzchen)

#################
## Pkgs needed ##
#################
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(rlang))
# conflicted::conflict_prefer("filter", "dplyr")

#' filterByDomains
#'
#' @author Samuel Chen, Janani Ravi
#' @description filterByDomains filters a data frame by identifying exact domain matches
#' and either keeping or removing rows with the identified domain
#'
#' @param prot Dataframe to filter
#' @param column Column to search for domains in (DomArch column)
#' @param doms_keep Vector of domains that must be identified within column in order for
#' observation to be kept
#' @param doms_remove Vector of domains that, if found within an observation, will be removed
#' @param ignore.case Should the matching be non case sensitive
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_replace_all
#' @importFrom rlang abort sym
#'
#' @return Filtered data frame
#' @note There is no need to make the domains 'regex safe', that will be handled by this function
#' @export
#'
#' @examples
#' \dontrun{
#' filterByDomains()
#' }
filterByDomains <- function(prot, column = "DomArch", doms_keep = c(), doms_remove = c(),
    ignore.case = FALSE) {
    # Only rows with a domain in doms_keep will be kept
    # Any row containing a domain in doms_remove will be removed

    # ^word$|(?<=\+)word$|(?<=\+)word(?=\+)|word(?=\+)
    
    # Check if prot is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        abort(paste("Error: The specified column '", column, "' does not exist 
                   in the data frame.", sep = ""))
    }
    
    # If doms_keep or doms_remove are not provided, inform the user
    if (length(doms_keep) == 0 && length(doms_remove) == 0) {
        warning("Warning: No domains specified to keep or remove. Returning the
                original data frame.")
    }

    # Make regex safe
    doms_keep <- str_replace_all(string = doms_keep, pattern = "\\(", replacement = "\\\\(")
    doms_keep <- str_replace_all(string = doms_keep, pattern = "\\)", replacement = "\\\\)")
    doms_keep <- str_replace_all(string = doms_keep, pattern = "\\+", replacement = "\\\\+")
    doms_keep <- str_replace_all(string = doms_keep, pattern = "\\_", replacement = "\\\\_")
    doms_keep <- str_replace_all(string = doms_keep, pattern = "\\?", replacement = "\\\\?")

    doms_remove <- str_replace_all(string = doms_remove, pattern = "\\(", replacement = "\\\\(")
    doms_remove <- str_replace_all(string = doms_remove, pattern = "\\)", replacement = "\\\\)")
    doms_remove <- str_replace_all(string = doms_remove, pattern = "\\+", replacement = "\\\\+")
    doms_remove <- str_replace_all(string = doms_remove, pattern = "\\_", replacement = "\\\\_")
    doms_remove <- str_replace_all(string = doms_remove, pattern = "\\?", replacement = "\\\\?")

    col <- sym(column)

    if (length(doms_keep) != 0) {
        keep_regex <- paste0(
            "^", doms_keep,
            "$|(?<=\\+)", doms_keep,
            "$|(?<=\\+)", doms_keep,
            "(?=\\+)|^", doms_keep,
            "(?=\\+)"
        ) %>%
            paste0(collapse = "|")
        prot <- prot %>% filter(grepl(keep_regex, {{ col }}, ignore.case = ignore.case, perl = T))
    }

    if (length(doms_remove) != 0) {
        remove_regex <- paste0(
            "^", doms_remove,
            "$|(?<=\\+)", doms_remove,
            "$|(?<=\\+)", doms_remove,
            "(?=\\+)|^", doms_remove,
            "(?=\\+)"
        ) %>%
            paste0(collapse = "|")
        prot <- prot %>% filter(!grepl(remove_regex, {{ col }}, ignore.case = ignore.case, perl = T))
    }

    return(prot)
}

###########################
## COUNTS of DAs and GCs ##
## Before/after break up ##
###########################

#' countByColumn
#' @description
#' Function to obtain element counts (DA, GC)
#' 
#' @param prot A data frame containing the dataset to analyze, typically with 
#' multiple columns including the one specified by the `column` parameter.
#' @param column A character string specifying the name of the column to analyze. 
#' The default is "DomArch".
#' @param min.freq An integer specifying the minimum frequency an element must 
#' have to be included in the output. Default is 1.
#'
#' @importFrom dplyr arrange as_tibble filter select
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{`column`}{The unique elements from the specified column 
#'   (e.g., "DomArch").}
#'   \item{`freq`}{The frequency of each element, i.e., the number of times 
#'   each element appears in the specified column.}
#' }
#' The tibble is filtered to only include elements that have a frequency 
#' greater than or equal to `min.freq` and does not include elements with `NA` 
#' values or those starting with a hyphen ("-").
#' @export
#'
#' @examples
#' \dontrun{
#' countByColumn(prot = my_data, column = "DomArch", min.freq = 10)
#' }
countByColumn <- function(prot = prot, column = "DomArch", min.freq = 1) {
    
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        abort(paste("Error: The specified column '", column, "' does not exist in
                   the data frame.", sep = ""))
    }
    
    # Check if min.freq is a positive integer
    if (!is.numeric(min.freq) || length(min.freq) != 1 || min.freq < 1 || 
        floor(min.freq) != min.freq) {
        abort("Error: 'min.freq' must be a positive integer.")
    }
    counts <- prot %>%
        select(column) %>%
        table() %>%
        as_tibble() %>%
        `colnames<-`(c(column, "freq")) %>%
        filter(!grepl("^-$", column)) %>% # remove "-"
        filter(!is.na(column)) %>%
        arrange(-freq) %>%
        filter(freq >= min.freq)
    return(counts)
}

#' elements2Words
#'
#' @description
#' Break string ELEMENTS into WORDS for domain architecture (DA) and genomic
#' context (GC)
#'
#' @param prot A dataframe containing the dataset to analyze. The specified 
#' `column` contains the string elements to be processed.
#' @param column A character string specifying the name of the column to analyze. 
#' Default is "DomArch".
#' @param conversion_type A character string specifying the type of conversion. 
#' Two options are available:
#' \describe{
#'   \item{`da2doms`}{Convert domain architectures into individual domains by 
#'   replacing `+` symbols with spaces.}
#'   \item{`gc2da`}{Convert genomic context into domain architectures by
#'    replacing directional symbols (`<-`, `->`, and `|`) with spaces.}
#' }
#'
#' @importFrom dplyr pull
#' @importFrom stringr str_replace_all
#'
#' @return A single string where elements are delimited by spaces. The function 
#' performs necessary substitutions based on the `conversion_type` and cleans up 
#' extraneous characters like newlines, tabs, and multiple spaces.
#'
#' @examples
#' \dontrun{
#' tibble::tibble(DomArch = c("aaa+bbb", 
#' "a+b", "b+c", "b-c")) |> elements2Words()
#' }
#'
elements2Words <- function(prot, column = "DomArch", conversion_type = "da2doms") {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        abort(paste("Error: The specified column '", column, "' does not exist in 
                   the data frame.", sep = ""))
    }
    
    # Check for valid conversion_type values
    valid_types <- c("da2doms", "doms2da")
    if (!conversion_type %in% valid_types) {
        abort(paste("Error: Invalid 'conversion_type'. Must be one of:", 
                   paste(valid_types, collapse = ", ")))
    }
    
    z1 <- prot %>%
        dplyr::pull(column) %>%
        str_replace_all("\\,", " ") %>%
        str_replace_all("\"", " ")
    switch(conversion_type,
        da2doms = {
            z2 <- z1 %>%
                str_replace_all("\\+", " ")
        },
        gc2da = {
            z2 <- z1 %>%
                str_replace_all("\\<-", " ") %>%
                str_replace_all("-\\>", " ") %>%
                str_replace_all("\\|", " ")
        }
    )
    # str_replace_all("^c\\($", " ") %>%		# remove "c("
    # str_replace_all("\\)$", " ") %>%			# remove ")"
    # str_replace_all("\\(s\\)"," ") %>%		# Ignoring repeats
    # str_replace_all("-"," ") %>%
    ## replace \n, \r, \t
    z3 <- z2 %>%
        str_replace_all("\n", " ") %>%
        str_replace_all("\r", " ") %>%
        str_replace_all("\t", " ") %>%
        # reduce spaces with length 2 or greater to a single space
        str_replace_all("\\s{2,}", " ")
    z3 <- z3 |> paste(collapse = " ")
    return(z3)
}

#' words2WordCounts
#'
#' @description
#' Get word counts (wc) [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
#'
#' @param string A character string containing the elements (words) to count. 
#' This would typically be a space-delimited string representing domain 
#' architectures or genomic contexts.
#'
#' @importFrom dplyr as_tibble filter arrange
#' @importFrom stringr str_replace_all
#'
#' @return A tibble (tbl_df) with two columns: 
#' \describe{
#'   \item{`words`}{A column containing the individual words 
#'   (domains or domain architectures).}
#'   \item{`freq`}{A column containing the frequency counts for each word.}
#' }
#' 
#'
#' @examples
#' \dontrun{
#' tibble::tibble(DomArch = c("aaa+bbb", "a+b", "b+c", "b-c")) |>
#'     elements2Words() |>
#'     words2WordCounts()
#' }
#'
words2WordCounts <- function(string) {
    # Check if 'string' is a character vector of length 1
    if (!is.character(string) || length(string) != 1) {
        abort("Error: 'string' must be a single character vector.")
    }
    
    df_word_count <- string %>%
        # reduce spaces with length 2 or greater to a single space
        str_replace_all("\\s{2,}", " ") %>%
        paste(collapse = " ") %>%
        strsplit(" ") %>%
        # filter(grepl(query.list[j], Query)) %>% # Create separate WCs for each Query
        # select(DA.wc) %>%
        table() %>%
        as_tibble() %>%
        `colnames<-`(c("words", "freq")) %>%
        ## filter out 'spurious-looking' domains
        filter(!grepl(" \\{n\\}", words)) %>%
        filter(!grepl("^c\\($", words)) %>% # remove "c("
        filter(!grepl("^\\)$", words)) %>% # remove ")"
        filter(!grepl("^-$", words)) %>% # remove "-"
        filter(!grepl("^$", words)) %>% # remove empty rows
        filter(!grepl("^\\?$", words)) %>% # remove "?"
        filter(!grepl("^\\?\\*$", words)) %>% # remove "?*"
        filter(!grepl("^tRNA$", words)) %>% # remove "tRNA"
        filter(!grepl("^ncRNA$", words)) %>% # remove "ncRNA"
        filter(!grepl("^rRNA$", words)) %>% # remove "rRNA"
        filter(!grepl("^X$|^X\\(s\\)$", words)) %>% # remove "X" and "X(s)"

        # filter(!grepl("\\*", words)) %>%			# Remove/Keep only Query
        arrange(-freq)
    return(df_word_count)
}

#' filterByFrequency
#' @description
#' Function to filter based on frequencies
#' 
#' @param x A tibble (tbl_df) containing at least two columns: one for 
#' elements (e.g., `words`) and one for their frequency (e.g., `freq`).
#' @param min.freq A numeric value specifying the minimum frequency threshold. 
#' Only elements with frequencies greater than or equal to this value will be 
#' retained.
#'
#' @return A tibble with the same structure as `x`, but filtered to include 
#' only rows where the frequency is greater than or equal to `min.freq`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' filterByFrequency()
#' }
filterByFrequency <- function(x, min.freq) {
    
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    
    # Check if 'min.freq' is a positive integer
    if (!is.numeric(min.freq) || length(min.freq) != 1 || min.freq < 1 || 
        floor(min.freq) != min.freq) {
        abort("Error: 'min.freq' must be a positive integer.")
    }
    
    # Check if the 'freq' column exists in the data frame
    if (!"freq" %in% names(x)) {
        abort("Error: The data frame must contain a 'freq' column.")
    }
    x %>%
        filter(freq >= min.freq)
}

#########################
## SUMMARY FUNCTIONS ####
#########################
#' MolEvolvR Summary
#' @name MolEvolvR_summary
#' @description
#' A collection of summary functions for the MolEvolvR package.
#' 
NULL

#' summarizeByLineage
#'
#' @param prot A dataframe or tibble containing the data.
#' @param column A string representing the column to be summarized 
#' (e.g., `DomArch`). Default is "DomArch".
#' @param by A string representing the grouping column (e.g., `Lineage`). 
#' Default is "Lineage".
#' @param query A string specifying the query pattern for filtering the target 
#' column. Use "all" to skip filtering and include all rows.
#'
#' @importFrom dplyr arrange filter group_by summarise
#' @importFrom rlang sym
#'
#' @return A tibble summarizing the counts of occurrences of elements in 
#' the `column`, grouped by the `by` column. The result includes the number 
#' of occurrences (`count`) and is arranged in descending order of count.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' tibble(DomArch = c("a+b", "a+b", "b+c", "a+b"), Lineage = c("l1", "l1", "l1", "l2")) |>
#'     summarizeByLineage(query = "all")
#' }
#'
summarizeByLineage <- function(prot = "prot", column = "DomArch", by = "Lineage",
    query) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        abort(paste("Error: The specified column '", column, "' does not exist in 
                   the data frame.", sep = ""))
    }
    
    # Check if the 'by' column exists in the data frame
    if (!by %in% names(prot)) {
        abort(paste("Error: The specified 'by' column '", by, "' does not exist 
                   n the data frame.", sep = ""))
    }
    
    column <- sym(column)
    by <- sym(by)
    if (query == "all") {
        prot <- prot
    } else {
        prot <- prot %>% filter(grepl(
            pattern = query, x = {{ column }},
            ignore.case = T
        ))
    }
    prot %>%
        filter(!grepl("^-$", {{ column }})) %>%
        group_by({{ column }}, {{ by }}) %>%
        summarise(count = n()) %>% # , bin=as.numeric(as.logical(n()))
        arrange(desc(count))
}


#' summarizeDomArch_ByLineage
#'
#' @description
#' Function to summarize and retrieve counts by Domains & Domains+Lineage
#'
#'
#' @param x A dataframe or tibble containing the data. It must have columns 
#' named `DomArch` and `Lineage`.
#'
#' @importFrom dplyr arrange count desc filter group_by summarise
#' @importFrom rlang .data
#'
#' @return A tibble summarizing the counts of unique domain architectures 
#' (`DomArch`) per lineage (`Lineage`). The resulting table contains three 
#' columns: `DomArch`, `Lineage`, and `count`, which indicates the frequency 
#' of each domain architecture for each lineage. The results are arranged in 
#' descending order of `count`.
#' @rdname MolEvolvR_summary
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summarizeDomArch_ByLineage(data1)
#' }
summarizeDomArch_ByLineage <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    
    # Check if required columns exist in the data frame
    required_columns <- c("DomArch", "Lineage")
    missing_columns <- setdiff(required_columns, names(x))
    
    if (length(missing_columns) > 0) {
        abort(paste("Error: The following required columns are 
                   missing:", paste(missing_columns, collapse = ", ")))
    }
    x %>%
        filter(!grepl("^-$", .data$DomArch)) %>%
        group_by(.data$DomArch, .data$Lineage) %>%
        summarise(count = n()) %>% # , bin=as.numeric(as.logical(n()))
        arrange(desc(count))
}


#' summarizeDomArch
#'
#' @description
#' Function to retrieve counts of how many lineages a DomArch appears in
#'
#' @param x A dataframe or tibble containing the data. It must have a column 
#' named `DomArch` and a count column, such as `count`, which represents the 
#' occurrences of each architecture in various lineages.
#'
#' @importFrom dplyr arrange group_by filter summarise desc
#' @importFrom rlang .data
#'
#' @return A tibble summarizing each unique `DomArch`, along with the following 
#' columns:
#' - `totalcount`: The total occurrences of each `DomArch` across all lineages.
#' - `totallin`: The total number of unique lineages in which each `DomArch` 
#' appears.
#' The results are arranged in descending order of `totallin` and `totalcount`.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @examples
#' \dontrun{
#' summarizeDomArch(data1)
#' }
summarizeDomArch <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    x %>%
        group_by(.data$DomArch) %>%
        summarise(
            totalcount = sum(.data$count), 
            totallin = n()
        ) %>%
        arrange(desc(.data$totallin), desc(.data$totalcount)) %>%
        filter(!grepl(" \\{n\\}", .data$DomArch)) %>%
        filter(!grepl("^-$", .data$DomArch))
}

#' summarizeGenContext_ByDomArchLineage
#'
#' @param x A dataframe or tibble containing the data. It must have columns 
#' named `GenContext`, `DomArch`, and `Lineage`.
#'
#' @importFrom dplyr arrange desc filter group_by n summarise
#' @importFrom rlang .data
#'
#' @return A tibble summarizing each unique combination of `GenContext`, 
#' `DomArch`, and `Lineage`, along with the following columns:
#' - `GenContext`: The genomic context for each entry.
#' - `DomArch`: The domain architecture for each entry.
#' - `Lineage`: The lineage associated with each entry.
#' - `count`: The total number of occurrences for each combination of 
#' `GenContext`, `DomArch`, and `Lineage`.
#'
#' The results are arranged in descending order of `count`.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @examples
#' \dontrun{
#' summarizeGenContext_ByDomArchLineage(your_data)
#' }
summarizeGenContext_ByDomArchLineage <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    x %>%
        filter(!grepl("^-$", .data$GenContext)) %>%
        filter(!grepl("^-$", .data$DomArch)) %>%
        filter(!grepl("^-$", .data$Lineage)) %>%
        filter(!grepl("^NA$", .data$DomArch)) %>%
        group_by(.data$GenContext, .data$DomArch, .data$Lineage) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
}

#' summarizeGenContext_ByLineage
#'
#' @param x A dataframe or tibble containing the data. It must have columns 
#' named `GenContext`, `DomArch`, and `Lineage`.
#'
#' @importFrom dplyr arrange desc filter group_by n summarise
#' @importFrom rlang .data
#'
#' @return A tibble summarizing each unique combination of `GenContext` and `Lineage`, 
#' along with the count of occurrences. The results are arranged in descending order of count.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @examples
#' \dontrun{
#' summarizeGenContext_ByLineage(your_data)
#' }
summarizeGenContext_ByLineage <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    x %>%
        filter(!grepl("^-$", .data$GenContext)) %>%
        filter(!grepl("^-$", .data$DomArch)) %>%
        filter(!grepl("^-$", .data$Lineage)) %>%
        filter(!grepl("^NA$", .data$DomArch)) %>%
        group_by(.data$GenContext, .data$Lineage) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
}

#' summarizeGenContext
#'
#' @param x A dataframe or tibble containing the data. It must have columns 
#' named `GenContext`, `DomArch`, `Lineage`, and `count`.
#'
#' @importFrom dplyr arrange desc filter group_by n_distinct summarise
#' @importFrom rlang .data
#'
#' @return A tibble summarizing each unique `GenContext`, along with the following columns:
#' - `totalcount`: The total count for each `GenContext`.
#' - `totalDA`: The number of distinct `DomArch` for each `GenContext`.
#' - `totallin`: The number of distinct `Lineage` for each `GenContext`.
#'
#' The results are arranged in descending order of `totalcount`, `totalDA`, and `totallin`.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @examples
#' \dontrun{
#' summarizeGenContext(data1)
#' }
summarizeGenContext <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        abort("Error: 'x' must be a data frame.")
    }
    x %>%
        group_by(.data$GenContext) %>%
        summarise(
            totalcount = sum(.data$count),
            totalDA = n_distinct(.data$DomArch),
            totallin = n_distinct(.data$Lineage)
        ) %>%
        arrange(desc(.data$totalcount), desc(.data$totalDA), desc(.data$totallin)) %>%
        filter(!grepl(" \\{n\\}", .data$GenContext)) %>%
        filter(!grepl("^-$", .data$GenContext))
}


##################
#' totalGenContextOrDomArchCounts
#'
#' @description
#' Creates a data frame with a totalcount column
#'
#' This function is designed to sum the counts column by either Genomic Context or Domain Architecture and creates a totalcount column from those sums.
#'
#'
#' @param prot  A data frame that must contain columns:
#' \itemize{\item Either 'GenContext' or 'DomArch.norep' \item count}
#' @param column Character. The column to summarize, default is "DomArch".
#' @param lineage_col Character. The name of the lineage column, default is 
#' "Lineage".
#' @param cutoff Numeric. Cutoff for total count. Counts below this cutoff value 
#' will not be shown. Default is 0.
#' @param RowsCutoff Logical. If TRUE, filters based on cumulative percentage 
#' cutoff. Default is FALSE.
#' @param digits Numeric. Number of decimal places for percentage columns. 
#' Default is 2.
#'
#'
#' @importFrom dplyr arrange distinct filter group_by left_join mutate select summarise ungroup
#' @importFrom rlang as_string sym
#'
#' @return A data frame with the following columns:
#' - `{{ column }}`: Unique values from the specified column.
#' - `totalcount`: The total count of occurrences for each unique value in 
#' the specified column.
#' - `IndividualCountPercent`: The percentage of each `totalcount` relative to 
#' the overall count.
#' - `CumulativePercent`: The cumulative percentage of total counts.
#' @rdname MolEvolvR_summary
#' @export
#'
#' @note Please refer to the source code if you have alternate file formats and/or
#' column names.
#'
#' @examples
#' \dontrun{
#' totalGenContextOrDomArchCounts(pspa - gc_lin_counts, 0, "GC")
#' }
totalGenContextOrDomArchCounts <- function(prot, column = "DomArch", lineage_col = "Lineage",
    cutoff = 90, RowsCutoff = FALSE, digits = 2
    # type = "GC"
) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified columns exist in the data frame
    required_columns <- c(column, lineage_col)
    missing_columns <- setdiff(required_columns, names(prot))
    
    if (length(missing_columns) > 0) {
        abort(paste("Error: The following required columns are missing:", 
                   paste(missing_columns, collapse = ", ")))
    }
    
    # Check that cutoff is a numeric value between 0 and 100
    if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff < 0 || cutoff > 100) {
        abort("Error: 'cutoff' must be a numeric value between 0 and 100.")
    }
    
    # Check that digits is a non-negative integer
    if (!is.numeric(digits) || length(digits) != 1 || digits < 0 || 
        floor(digits) != digits) {
        abort("Error: 'digits' must be a non-negative integer.")
    }
    
    column <- sym(column)

    prot <- select(prot, {{ column }}, {{ lineage_col }}) %>%
        filter(!is.na({{ column }}) & !is.na({{ lineage_col }})) %>%
        filter({{ column }} != "")

    prot <- summarizeByLineage(prot, column, by = lineage_col, query = "all")
    col_count <- prot %>%
        group_by({{ column }}) %>%
        summarise(totalcount = sum(count))

    total <- left_join(prot, col_count, by = as_string(column))

    sum_count <- sum(total$count)
    total <- total %>%
        mutate("IndividualCountPercent" = totalcount / sum_count * 100) %>%
        arrange(-totalcount, -count)

    cumm_percent <- total %>%
        select({{ column }}, totalcount) %>%
        distinct() %>%
        mutate("CumulativePercent" = 0)
    total_counter <- 0
    for (x in length(cumm_percent$totalcount):1) {
        total_counter <- total_counter + cumm_percent$totalcount[x]
        cumm_percent$CumulativePercent[x] <- total_counter / sum_count * 100
    }

    cumm_percent <- cumm_percent %>% select(CumulativePercent, {{ column }})

    total <- total %>% left_join(cumm_percent, by = as_string(column))

    # Round the percentage columns
    total$CumulativePercent <- total$CumulativePercent %>% round(digits = digits)
    total$IndividualCountPercent <- total$IndividualCountPercent %>% round(digits = digits)

    if (RowsCutoff) {
        # If total counts is being used for plotting based on number of rows,
        # don't include other observations that fall below the cummulative percent cutoff
        # , but that have the same 'totalcount' number as the cutoff observation
        total <- total %>% filter(CumulativePercent >= 100 - cutoff - .0001)
        return(total)
    }

    # Include observations that fall below the cummulative percent cutoff,
    # but that have the same 'totalcount' as the cutoff observation
    t <- total %>% filter(CumulativePercent >= 100 - cutoff)
    if (length(t) == 0) {
        cutoff_count <- 0
    } else {
        cutoff_count <- t$totalcount[nrow(t)]
    }

    total <- total %>%
        filter(totalcount >= cutoff_count) %>%
        ungroup()

    return(total)
}



# total_counts_by_query <- function(query_data, queries, colname,cutoff, RemoveAstrk = F)
# {
#   ## Get the total counts by the Queries.
#
#   lineage_by_query <- function(data, query, column, by){
#     # Function to filter data by a query and summarise it by lineage
#
#     column <- sym(column); by <- sym(by)
#
#     # filter the protein by the query
#     data <- data %>% filter(grepl(pattern=query, x={{column}},
#                                   ignore.case=T)) %>% select({{by}})
#
#     data$Query <- query
#
#     data <- data %>% filter(!grepl("^-$", {{by}})) %>%
#       group_by(Query, {{by}}) %>%
#       summarise(count=n()) %>%
#       arrange(desc(count))
#     return(data)
#   }
#
#
#   col <- sym(colname)
#
#   # query_data contains all rows that possess a lineage
#   query_data <- query_data %>% filter(grepl("a", Lineage))
#
#   query_data <- shorten_lineage(query_data, "Lineage")
#
#   if(RemoveAstrk){
#     # Remove Asterisk (*) if necessary
#     query_data <- query_data %>% remove_astrk(colname = colname)
#   }
#
#   query_lin_counts = data.frame("Query" = character(0), "Lineage" = character(0), "count"= integer())
#
#   for(q in queries){
#     # iterate over queries and do filtering/Lineage by query
#     query_lin <- lineage_by_query(data = query_data, query = q, column = {{col}}, by = "Lineage")
#     query_lin_counts <- dplyr::union(query_lin_counts, query_lin)
#   }
#
#   # Total counts of each lineage
#   lin_count_totals <- query_lin_counts %>%group_by(Lineage) %>% summarize(total_c = sum(count)) %>% arrange(total_c)
#   sum_lin <- sum(lin_count_totals$total_c)
#
#   # Calculate the percent each lineage makes up
#   lin_count_totals$IndividualPercent = lin_count_totals$total_c/sum_lin * 100
#   lin_count_totals <- lin_count_totals %>% mutate("CumulativePercent"=0)
#   total_counter = 0
#   for(x in 1:length(lin_count_totals$IndividualPercent)){
#     total_counter = total_counter + lin_count_totals$IndividualPercent[x]
#     lin_count_totals$CumulativePercent[x] = total_counter
#   }
#
#   query_total_counts <- left_join(query_lin_counts, lin_count_totals, by = "Lineage")
#
#   query_total_counts <- query_total_counts%>% group_by(total_c,count) %>% arrange(-total_c, -count)
#   # Get lineages that are above the cutoff percentage value
#   count_cutoff <- (query_total_counts %>% filter(CumulativePercent >= (100-cutoff)))
#
#   query_total_counts <- (query_total_counts %>% filter(total_c >= count_cutoff$total_c[nrow(count_cutoff)]))
#
#   # Round to 3 digits
#   query_total_counts$CumulativePercent <- query_total_counts$CumulativePercent %>% round(digits = 3)
#   query_total_counts$IndividualPercent <- query_total_counts$IndividualPercent %>% round(digits = 3)
#
#   return(query_total_counts)
# }






#' findParalogs
#'
#' @description
#' Creates a data frame of paralogs.
#'
#' @param prot A data frame filtered by a Query, containing columns Species and Lineage
#'
#' @importFrom dplyr arrange count distinct filter group_by right_join select
#' @importFrom purrr map
#'
#' @return returns a dataframe containing paralogs and the counts.
#' @export
#'
#' @note Please refer to the source code if you have alternate file formats and/or
#' column names.
#'
#' @examples
#' \dontrun{
#' findParalogs(pspa)
#' }
findParalogs <- function(prot) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        abort("Error: 'prot' must be a data frame.")
    }
    
    # Remove eukaryotes
    prot <- prot %>% filter(!grepl("^eukaryota", Lineage))
    paralogTable <- prot %>%
        group_by(GCA_ID) %>%
        count(DomArch) %>%
        filter(n > 1 & GCA_ID != "-") %>%
        arrange(-n) %>%
        distinct()
    # %>% ccNums" = filter(prot, grepl(GCA_ID, prot))$AccNum)
    paralogTable$AccNums <- map(paralogTable$GCA_ID, function(x) filter(prot, grepl(x, GCA_ID))$AccNum)
    colnames(paralogTable)[colnames(paralogTable) == "n"] <- "Count"
    ### Merge with columns: AccNum,TaxID, and GCA/ Species?
    paralogTable <- prot %>%
        select(Species, GCA_ID, Lineage) %>%
        right_join(paralogTable, by = c("GCA_ID")) %>%
        filter(!is.na(Count)) %>%
        arrange(-Count) %>%
        select(-GCA_ID) %>%
        distinct()
    return(paralogTable)
}


##################################
## Descriptions for functions ####
##################################
# ## counts: Function to obtain element counts (DA, GC)
# cat("Counts DAs and GCs.
# \nFor e.g.:
# query.sub$DomArch.norep %>%
#   counts(n)
# query.sub$GenContext %>%
# counts(n)")

# ## elements2Words: Function to break up ELEMENTS to WORDS for DA and GC
# cat("Converting DA to domains and GC to DAs.\n2 switches: da2doms and gc2da
# \nFor e.g.:
# query.sub$DA.doms <- query.sub$DomArch.norep %>%
#   elements2Words(\"da2doms\")
# query.sub$GC.da <- query.sub$GenContext %>%
# 	elements2Words(\"gc2da\")")


# ## words2WordCounts: Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
# cat("Word counts for broken up domains from DAs and DAs from GCs.
# \nFor e.g.:
# DA.doms.wc <- query.sub$DA.doms %>%
#   words2WordCounts()")