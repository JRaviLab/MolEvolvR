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

#' Filter by Domains
#'
#' @author Samuel Chen, Janani Ravi
#' @description filter_by_doms filters a data frame by identifying exact domain matches
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
#' @importFrom rlang sym
#'
#' @return Filtered data frame
#' @note There is no need to make the domains 'regex safe', that will be handled by this function
#' @export
#'
#' @examples
#' \dontrun{
#' filter_by_doms()
#' }
filter_by_doms <- function(prot, column = "DomArch", doms_keep = c(), doms_remove = c(),
    ignore.case = FALSE) {
    # Only rows with a domain in doms_keep will be kept
    # Any row containing a domain in doms_remove will be removed

    # ^word$|(?<=\+)word$|(?<=\+)word(?=\+)|word(?=\+)
    
    # Check if prot is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        stop(paste("Error: The specified column '", column, "' does not exist 
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
## Function to obtain element counts (DA, GC)
#' Count Bycol
#'
#' @param prot
#' @param column
#' @param min.freq
#'
#' @importFrom dplyr arrange as_tibble filter select
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' count_bycol()
#' }
count_bycol <- function(prot = prot, column = "DomArch", min.freq = 1) {
    
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        stop(paste("Error: The specified column '", column, "' does not exist in
                   the data frame.", sep = ""))
    }
    
    # Check if min.freq is a positive integer
    if (!is.numeric(min.freq) || length(min.freq) != 1 || min.freq < 1 || 
        floor(min.freq) != min.freq) {
        stop("Error: 'min.freq' must be a positive integer.")
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

#' Elements 2 Words
#'
#' @description
#' Break string ELEMENTS into WORDS for domain architecture (DA) and genomic
#' context (GC)
#'
#' @param prot [dataframe]
#' @param column [string] column name
#' @param conversion_type [string] type of conversion: 'da2doms': domain architectures to
#' domains. 'gc2da' genomic context to domain architectures
#'
#' @importFrom dplyr pull
#' @importFrom stringr str_replace_all
#'
#' @return [string] with words delimited by spaces
#'
#' @examples
#' \dontrun{
#' tibble::tibble(DomArch = c("aaa+bbb", "a+b", "b+c", "b-c")) |> elements2words()
#' }
#'
elements2words <- function(prot, column = "DomArch", conversion_type = "da2doms") {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        stop(paste("Error: The specified column '", column, "' does not exist in 
                   the data frame.", sep = ""))
    }
    
    # Check for valid conversion_type values
    valid_types <- c("da2doms", "doms2da")
    if (!conversion_type %in% valid_types) {
        stop(paste("Error: Invalid 'conversion_type'. Must be one of:", 
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

#' Words 2 Word Counts
#'
#' @description
#' Get word counts (wc) [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
#'
#' @param string
#'
#' @importFrom dplyr as_tibble filter
#'
#' @return [tbl_df] table with 2 columns: 1) words & 2) counts/frequency
#'
#' @examples
#' \dontrun{
#' tibble::tibble(DomArch = c("aaa+bbb", "a+b", "b+c", "b-c")) |>
#'     elements2words() |>
#'     words2wc()
#' }
#'
words2wc <- function(string) {
    # Check if 'string' is a character vector of length 1
    if (!is.character(string) || length(string) != 1) {
        stop("Error: 'string' must be a single character vector.")
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
## Function to filter based on frequencies
#' Filter Frequency
#'
#' @param x
#' @param min.freq
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' filter_freq()
#' }
filter_freq <- function(x, min.freq) {
    
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    
    # Check if 'min.freq' is a positive integer
    if (!is.numeric(min.freq) || length(min.freq) != 1 || min.freq < 1 || 
        floor(min.freq) != min.freq) {
        stop("Error: 'min.freq' must be a positive integer.")
    }
    
    # Check if the 'freq' column exists in the data frame
    if (!"freq" %in% names(x)) {
        stop("Error: The data frame must contain a 'freq' column.")
    }
    x %>%
        filter(freq >= min.freq)
}

#########################
## SUMMARY FUNCTIONS ####
#########################
#' Summarize by Lineage
#'
#' @param prot
#' @param column
#' @param by
#' @param query
#'
#' @importFrom dplyr arrange filter group_by summarise
#' @importFrom rlang sym
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' tibble(DomArch = c("a+b", "a+b", "b+c", "a+b"), Lineage = c("l1", "l1", "l1", "l2")) |>
#'     summarize_bylin(query = "all")
#' }
#'
summarize_bylin <- function(prot = "prot", column = "DomArch", by = "Lineage",
    query) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified column exists in the data frame
    if (!column %in% names(prot)) {
        stop(paste("Error: The specified column '", column, "' does not exist in 
                   the data frame.", sep = ""))
    }
    
    # Check if the 'by' column exists in the data frame
    if (!by %in% names(prot)) {
        stop(paste("Error: The specified 'by' column '", by, "' does not exist 
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


#' summ.DA.byLin
#'
#' @description
#' Function to summarize and retrieve counts by Domains & Domains+Lineage
#'
#'
#' @param x
#'
#' @importFrom dplyr arrange count desc filter group_by summarise
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' summ.DA.byLin()
#' }
summ.DA.byLin <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    
    # Check if required columns exist in the data frame
    required_columns <- c("DomArch", "Lineage")
    missing_columns <- setdiff(required_columns, names(x))
    
    if (length(missing_columns) > 0) {
        stop(paste("Error: The following required columns are 
                   missing:", paste(missing_columns, collapse = ", ")))
    }
    ## Note: it is better to reserve dots for S3 Objects. Consider replacing '.' with '_'
    x %>%
        filter(!grepl("^-$", DomArch)) %>%
        group_by(DomArch, Lineage) %>%
        summarise(count = n()) %>% # , bin=as.numeric(as.logical(n()))
        arrange(desc(count))
}

## Function to retrieve counts of how many lineages a DomArch appears in
#' summ.DA
#'
#' @description
#' Function to retrieve counts of how many lineages a DomArch appears in
#'
#' @param x
#'
#' @importFrom dplyr arrange group_by filter summarise
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' summ.DA()
#' }
summ.DA <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    ## Note: it is better to reserve dots for S3 Objects. Consider replacing '.' with '_'
    x %>%
        group_by(DomArch) %>%
        summarise(totalcount = sum(count), totallin = n()) %>% # totallin=n_distinct(Lineage),
        arrange(desc(totallin), desc(totalcount)) %>%
        filter(!grepl(" \\{n\\}", DomArch)) %>%
        filter(!grepl("^-$", DomArch))
}

#' summ.GC.byDALin
#'
#' @param x
#'
#' @importFrom dplyr arrange desc filter group_by n summarise
#'
#' @return Define return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' summ.GC.byDALin
#' }
summ.GC.byDALin <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    ## Note: it is better to reserve dots for S3 Objects. Consider replacing '.' with '_'
    x %>%
        filter(!grepl("^-$", GenContext)) %>%
        filter(!grepl("^-$", DomArch)) %>%
        filter(!grepl("^-$", Lineage)) %>%
        filter(!grepl("^NA$", DomArch)) %>%
        group_by(GenContext, DomArch, Lineage) %>%
        summarise(count = n()) %>% # , bin=as.numeric(as.logical(n()))
        arrange(desc(count))
}

#' summ.GC.byLin
#'
#' @param x
#'
#' @importFrom dplyr arrange desc filter group_by n summarise
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' summ.GC.byLin()
#' }
summ.GC.byLin <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    ## Note: it is better to reserve dots for S3 Objects. Consider replacing '.' with '_'
    x %>%
        filter(!grepl("^-$", GenContext)) %>%
        filter(!grepl("^-$", DomArch)) %>%
        filter(!grepl("^-$", Lineage)) %>%
        filter(!grepl("^NA$", DomArch)) %>%
        group_by(GenContext, Lineage) %>% # DomArch.norep,
        summarise(count = n()) %>% # , bin=as.numeric(as.logical(n()))
        arrange(desc(count))
}

#' summ.GC
#'
#' @param x
#'
#' @importFrom dplyr arrange desc filter group_by n_distinct summarise
#'
#' @return Describe return, in detail
#' @export
#'
#' @examples
#' \dontrun{
#' summ.GC()
#' }
summ.GC <- function(x) {
    # Check if 'x' is a data frame
    if (!is.data.frame(x)) {
        stop("Error: 'x' must be a data frame.")
    }
    ## Note: it is better to reserve dots for S3 Objects. Consider replacing '.' with '_'
    x %>%
        group_by(GenContext) %>%
        summarise(
            totalcount = sum(count),
            totalDA = n_distinct(DomArch),
            totallin = n_distinct(Lineage)
        ) %>% # totallin=n_distinct(Lineage),
        arrange(desc(totalcount), desc(totalDA), desc(totallin)) %>%
        filter(!grepl(" \\{n\\}", GenContext)) %>%
        filter(!grepl("^-$", GenContext))
}


##################
#' Total Counts
#'
#' @description
#' Creates a data frame with a totalcount column
#'
#' This function is designed to sum the counts column by either Genomic Context or Domain Architecture and creates a totalcount column from those sums.
#'
#'
#' @param prot  A data frame that must contain columns:
#' \itemize{\item Either 'GenContext' or 'DomArch.norep' \item count}
#' @param column Character. The column to summarize
#' @param lineage_col
#' @param cutoff Numeric. Cutoff for total count. Counts below cutoff value will not be shown. Default is 0.
#' @param RowsCutoff
#' @param digits
#'
#' @importFrom dplyr arrange distinct filter group_by left_join mutate select summarise ungroup
#' @importFrom rlang as_string sym
#'
#' @return Define return, in detail
#' @export
#'
#' @note Please refer to the source code if you have alternate file formats and/or
#' column names.
#'
#' @examples
#' \dontrun{
#' total_counts(pspa - gc_lin_counts, 0, "GC")
#' }
total_counts <- function(prot, column = "DomArch", lineage_col = "Lineage",
    cutoff = 90, RowsCutoff = FALSE, digits = 2
    # type = "GC"
) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
    }
    
    # Check if the specified columns exist in the data frame
    required_columns <- c(column, lineage_col)
    missing_columns <- setdiff(required_columns, names(prot))
    
    if (length(missing_columns) > 0) {
        stop(paste("Error: The following required columns are missing:", 
                   paste(missing_columns, collapse = ", ")))
    }
    
    # Check that cutoff is a numeric value between 0 and 100
    if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff < 0 || cutoff > 100) {
        stop("Error: 'cutoff' must be a numeric value between 0 and 100.")
    }
    
    # Check that digits is a non-negative integer
    if (!is.numeric(digits) || length(digits) != 1 || digits < 0 || 
        floor(digits) != digits) {
        stop("Error: 'digits' must be a non-negative integer.")
    }
    
    column <- sym(column)

    prot <- select(prot, {{ column }}, {{ lineage_col }}) %>%
        filter(!is.na({{ column }}) & !is.na({{ lineage_col }})) %>%
        filter({{ column }} != "")

    prot <- summarize_bylin(prot, column, by = lineage_col, query = "all")
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






#' Find Paralogs
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
#' find_paralogs(pspa)
#' }
find_paralogs <- function(prot) {
    # Check if 'prot' is a data frame
    if (!is.data.frame(prot)) {
        stop("Error: 'prot' must be a data frame.")
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

# ## elements2words: Function to break up ELEMENTS to WORDS for DA and GC
# cat("Converting DA to domains and GC to DAs.\n2 switches: da2doms and gc2da
# \nFor e.g.:
# query.sub$DA.doms <- query.sub$DomArch.norep %>%
#   elements2words(\"da2doms\")
# query.sub$GC.da <- query.sub$GenContext %>%
# 	elements2words(\"gc2da\")")


# ## words2wc: Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
# cat("Word counts for broken up domains from DAs and DAs from GCs.
# \nFor e.g.:
# DA.doms.wc <- query.sub$DA.doms %>%
#   words2wc()")
