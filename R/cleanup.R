## Functions to clean up .op_ins_cls files by:
## Species, ClustName, DomArch, GenContext
## To create consistent names and take care of repeats & remove empty rows
## Created: Aug 11, 2017
## Modified: Dec 18, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
library(glue)
library(rlang) # needed?
conflicted::conflict_prefer("filter", "dplyr")

###########################
#### CLEANUP FUNCTIONS ####
###########################
remove_empty <- function(prot, by_column = "DomArch") {
  #' Remove empty rows by column
  #'
  #' Removes empty rows in the specified column.
  #'
  #' This function ...
  #' The original data frame is returned with the corresponding cleaned up column.
  #'
  #' @param prot A data frame containing 'DomArch', 'Species', 'GenContext', 'ClustName' columns.
  #' @param by_column Column by which empty rows should be removed to domain+domain -> domain(s).
  #' Default column is 'DomArch'. Can also take the following as input, 'Species', 'GenContext', 'ClustName'.
  #' @examples remove_empty(prot, "DomArch")

  # ?? Don't call other psp functions within these functions
  prot <- prot %>%
    as_tibble() %>%
    # filter(grepl("\\*", {{by_column}})) %>%		  # Keep only rows with Query (*) for GenContext
    filter(!grepl("^-$", {{ by_column }})) %>% # remove "-"
    filter(!grepl("^NA$", {{ by_column }})) %>% # remove "NA"
    filter(!grepl("^$", {{ by_column }})) # remove empty rows

  return(prot)
}

###########################
repeat2s <- function(prot, by_column = "DomArch", excluded_prots = c()) {
  #' Condense repeated domains
  #'
  #' Condenses repeated domains in the specified column.
  #'
  #' This function identifies repeated domains and condenses them to (s).
  #' ?? Certain domains can be removed according to an additional data frame.
  #' The original data frame is returned with the corresponding cleaned up column.
  #'
  #' @param prot A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.
  #' @param by_column Column in which repeats are condensed to domain+domain -> domain(s).
  #' Default column is 'DomArch'. Can also take the following as input, 'GenContext', 'ClustName'.
  #' @param excluded_prots Vector of strings that repeat2s should not reduce to (s). Defaults to c()
  #' @examples repeat2s(prot, "DomArch")


  # If there are strings that repeat2s should not affect, the pattern to search
  # for must be changed to exclude a search for those desired strings

  collapsed_prots <- paste0(excluded_prots, collapse = "\\s|")
  regex_exclude <- paste0("(?!", collapsed_prots, "\\s)")
  regex_identify_repeats <- paste0("(?i)", regex_exclude, "\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+")

  # !! FUNS is soft-deprecated. FIX!!!
  prot[, by_column] <- prot %>%
    pull(by_column) %>%
    str_replace_all(., pattern = "\\.", replacement = "_d_") %>%
    #  str_replace_all(., pattern = " ", replacement = "_s_") %>%
    str_replace_all(., pattern = " ", replacement = "_") %>%
    str_replace_all(.,
      pattern = "\\+",
      replacement = " "
    ) %>% # Use a different placeholder other than space
    str_replace_all(.,
      pattern = "-",
      replacement = "__"
    ) %>%
    str_replace_all(.,
      pattern = regex_identify_repeats,
      replacement = "\\1(s)"
    ) %>%
    str_replace_all(.,
      pattern = "__",
      replacement = "-"
    ) %>%
    str_replace_all(.,
      pattern = " ",
      replacement = "+"
    ) %>%
    # 			    str_replace_all(., pattern = "_s_", replacement = " ") %>%
    str_replace_all(., pattern = "_d_", replacement = ".")


  return(prot)
}


replaceQMs <- function(prot, by_column = "GenContext") {
  #' Replace consecutive '?' separated by '->', '<-' or '||' with 'X(s)'
  #' Replace '?' with 'X'
  #' @param prot DataTable to operate on
  #' @param by_column Column to operate on

  by <- sym(by_column)

  # Regex for finding repeated `?`
  regex_findRepeatedQ <- "(?i)(?!\\s)([?]+)(?:[\\|><\\s-]+\\1)+"

  # Replace Repeated ? with X(s) and ? with 'X'
  prot[, by_column] <- prot %>%
    pull({{ by }}) %>%
    str_replace_all(pattern = regex_findRepeatedQ, replacement = "X(s)") %>%
    str_replace_all(pattern = "\\?", replacement = "X")


  return(prot)
}


remove_astrk <- function(query_data, colname = "GenContext") {
  # Remove the asterisks from a column of data
  # Used for removing * from GenContext columns
  query_data[, colname] <- map(query_data[, colname], function(x) str_remove_all(x, pattern = "\\*"))

  return(query_data)
}


###########################
remove_tails <- function(prot, by_column = "DomArch",
                         keep_domains = FALSE) { # !! currently redundant
  #' Remove tails/singletons
  #'
  #'
  #' This function ...
  #' Certain low frequency domain architectures can be removed.
  #' The original data frame is returned with the corresponding cleaned up column.
  #'
  #' @param prot A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.
  #' @param by_column Default column is 'DomArch'. Can also take 'ClustName', 'GenContext' as input.
  #' @param keep_domains Default is False Keeps tail entries that contain the query domains.
  #' @examples remove_tails(prot, "DomArch")
  by_column <- sym(by_column)
  domain_count <- prot %>%
    group_by({{ by_column }}) %>%
    summarize(count = n())

  ## Identify tails
  tails <- domain_count %>% filter(count == 1)

  ## Domains_keep
  if (keep_domains) {
    # Keep tails with query domains
    # !! Insert line to read domains_keep

    # Contains all domains separated by "|"
    domains_for_grep <- paste(domains_keep$domains, collapse = "|")
    # Remove rows with no domains contained within domains_keep
    # Redundant for ClustName since we already set the filter to only these doms.
    tails <- tails %>%
      filter(!grepl(domains_for_grep, {{ by_column }}))
  }

  # Remove tails
  # tails <- tails %>% select({{by_column}}) %>% unlist()
  tails <- tails %>% pull({{ by_column }})
  prot <- prot %>% filter(!({{ by_column }} %in% tails))

  return(prot)
}

###########################
cleanup_species <- function(prot, remove_empty = FALSE) {
  #' Cleanup Species
  #'
  #' Cleans up the species column of a data frame by removing certain characters and rows.
  #'
  #' This function removes unneccessary characters from the 'Species' column.
  #' A cleaned up version of the data table is returned.
  #'
  #' @param prot A data frame that contains columns 'Species'.
  #' @param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'Species' are removed.
  #' Default is false.
  #' @examples cleanup_species(prot,TRUE)
  # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F

  # Create cleaned up Species column
  prot$Species <- prot$Species.orig %>%
    # remove dots after sp and str
    str_replace_all(coll("sp. ", TRUE), "sp ") %>%
    str_replace_all(coll("str. ", TRUE), "str ") %>%
    # remove special characters and brackets
    str_replace_all(coll(" = ", TRUE), " ") %>%
    str_replace_all(coll("-", TRUE), "") %>%
    str_replace_all(coll(".", TRUE), "") %>%
    str_replace_all(coll("(", TRUE), "") %>%
    str_replace_all(coll(")", TRUE), "") %>%
    str_replace_all(coll("[", TRUE), "") %>%
    str_replace_all(coll("]", TRUE), "") %>%
    str_replace_all(coll("\\", TRUE), "") %>%
    str_replace_all(coll("/", TRUE), "") %>%
    str_replace_all(coll("\'", TRUE), "") %>%
    # remove extra spaces
    str_replace_all(coll("  ", TRUE), " ")

  # !! CHECK !! Species vs Species_old
  if (remove_empty) {
    prot <- remove_empty(prot = prot, by_column = "Species")
  }

  return(prot)
}

######################
cleanup_clust <- function(prot,
                          domains_rename, domains_keep,
                          repeat2s = TRUE, remove_tails = FALSE,
                          remove_empty = FALSE) {
  #' Cleanup cluster file
  #'
  #' Cleans a cluster file by removing rows that do not contain the query in the cluster.
  #'
  #' This function removes irrelevant rows which do not contain the query protein within the ClustName column.
  #' The return value is the cleaned up data frame.
  #'
  #' @param prot A data frame that must contain columns Query and ClustName.
  #' @param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the corresponding replacement values in a column 'new'.
  #' @param domains_keep A data frame containing the domain names to be retained.
  #' @param repeat2s Boolean. If TRUE, repeated domains in 'ClustName' are condensed. Default is TRUE.
  #' @param remove_tails Boolean. If TRUE, 'ClustName' will be filtered based on domains to keep/remove. Default is FALSE.
  #' @param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'ClustName' are removed. Default is FALSE.

  #' @examples cleanup_clust(prot, TRUE, FALSE, domains_keep, domains_rename)

  # Create cleaned up ClustName column
  prot$ClustName <- prot$ClustName.orig

  ## Basic Cleanup
  # Remove '+' at the start and end, as well as consecuative '+'
  prot$ClustName <- gsub("\\++\\+", "\\+", prot$ClustName)
  prot$ClustName <- gsub("^\\+", "", prot$ClustName)
  prot$ClustName <- gsub("\\+$", "", prot$ClustName)

  ## Domains_rename
  for (x in 1:length(domains_rename$old)) {
    target <- domains_rename$old[x]
    replacement <- domains_rename$new[x]
    prot$ClustName <- prot$ClustName %>%
      str_replace_all(coll(target, TRUE), replacement)
  }

  ## Domains_keep
  # Character for greping for rows with domains_keep
  # Contains all domains separated by "|"
  domains_for_grep <- paste(domains_keep$domains, collapse = "|")
  # Remove rows with no domains contained within domains_keep
  prot <- prot %>% filter(grepl(domains_for_grep, ClustName))

  ## Optional parameters
  # Condense repeats
  if (repeat2s) {
    prot <- repeat2s(prot, by_column = "ClustName")
  }
  # Remove singletons
  # if(remove_tails){
  #  prot <- prot %>% filter(!grepl(".1$", ClustID))
  # }
  if (remove_tails) {
    prot <- remove_tails(prot, by_column = "ClustName")
  }
  # Remove empty rows
  if (remove_empty) {
    prot <- remove_empty(prot = prot, by_column = "ClustName")
  }


  # !!UNFIXED ISSUE!! Currently requires manual intervention!
  # SIG+TM+TM+... kind of architectures without explicit domain names are lost.
  # Need a way to take care of true hits that don't go by the expected domain name.
  return(prot)
}

##############################
cleanup_domarch <- function(prot, old = "DomArch.orig", new = "DomArch",
                            domains_keep, domains_rename,
                            repeat2s = TRUE, remove_tails = FALSE,
                            remove_empty = F,
                            domains_ignore = NULL) {
  #' Cleanup Domain Architectures
  #'
  #' Cleans the DomArch column by replacing/removing certain domains
  #'
  #' This function cleans the DomArch column of one data frame by renaming certain domains according to a second data frame.
  #' Certain domains can be removed according to an additional data frame.
  #' The original data frame is returned with the clean DomArchs column and the old domains in the DomArchs.old column.
  #'
  #' @param prot A data frame containing a 'DomArch' column
  #' @param domains_keep A data frame containing the domain names to be retained.
  #' @param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the
  #' corresponding replacement values in a column 'new'.
  #' @param repeat2s Boolean. If TRUE, repeated domains in 'DomArch' are condensed. Default is TRUE.
  #' @param remove_tails Boolean. If TRUE, 'ClustName' will be filtered based on domains to keep/remove. Default is FALSE.
  #' @param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'DomArch' are removed. Default is FALSE.
  #' @param domains_ignore A data frame containing the domain names to be removed in a column called 'domains'
  #' @examples cleanup_domarch(prot, TRUE, FALSE, domains_keep, domains_rename, domains_ignore=NULL)

  old_sym <- sym(old)
  new_sym <- sym(new)

  # Create cleaned up DomArch column
  # prot$DomArch <- prot$DomArch.orig
  prot[, new] <- prot[, old]

  ## Basic Cleanup
  # Remove '+' at the start and end, as well as consecuative '+'
  prot[, new] <- gsub("\\++\\+", "\\+", pull(prot, {{ new_sym }}))
  prot[, new] <- gsub("^\\+", "", pull(prot, {{ new_sym }}))
  prot[, new] <- gsub("\\+$", "", pull(prot, {{ new_sym }}))

  ## Domains_rename
  # Replace domains based on the domains_rename list

  if (!is.null(domains_rename) && length(domains_rename$old) != 0) {
    for (j in 1:length(domains_rename$old)) {
      prot[, new] <- str_replace_all(
        pull(prot, {{ new_sym }}),
        coll(as.vector(domains_rename$old[j]), TRUE),
        as.vector(domains_rename$new[j])
      )
    }
  }
  ## Domains_keep
  # Character for greping for rows with domains_keep
  # Contains all domains separated by "|"
  # domains_for_grep <- paste(domains_keep$domains, collapse = "|")
  # Remove rows with no domains contained within domains_keep
  # filter(grepl(domains_for_grep, DomArch))
  if (!is.null(domains_keep)) {
    prot <- prot %>% filter_by_doms(column = new, doms_keep = domains_keep$domains)
  }

  # ##!! NOT RUN !!
  # ## Domains_ignore
  # # Remove domains based on the domains_ignore list
  # # ?? Remove the domains or the rows? Check, please!
  # if( !is.null(domains_ignore)){
  #   for(j in 1:length(domains_ignore$domains)){
  #     prot$DomArch <- str_remove_all(prot$DomArch,
  #                                    as.vector(domains_ignore$domains[j]))
  #   }
  # }

  ## Optional parameters
  # Remove singletons
  if (remove_tails) {
    prot <- remove_tails(prot = prot, by_column = new)
  }
  # Condense repeats
  if (repeat2s) {
    ## Error in UseMethod("tbl_vars") : no applicable method for 'tbl_vars' applied to an object of class "character"
    prot <- repeat2s(prot = prot, by_column = new)
  }
  # Remove empty rows
  # ! FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  # ! Make a separate function of out of this?
  if (remove_empty) {
    prot <- remove_empty(prot = prot, by_column = new)
  }

  prot <- replaceQMs(prot, new)

  return(prot)
}

##########################
cleanup_gencontext <- function(prot, domains_rename = data.frame("old" = character(0), "new" = character(0), stringsAsFactors = F),
                               repeat2s = TRUE, remove_asterisk = TRUE) {
  #' Cleanup Genomic Contexts
  #'
  #' Cleans up the GenContext column of a data frame by removing certain characters and rows.
  #'
  #' This function removes empty rows based on the 'GenContext' column.
  #' A cleaned up version of the data table is returned.
  #'
  #' @param prot A data frame that contains columns 'GenContext.orig'
  #' @param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the replacement in a column 'new'.
  #' Defaults to an empty data frame with a new and old column such that non of the domains will be renamed
  #' @param repeat2s Boolean. If TRUE, repeated domains in 'GenContext' are condensed. Default is TRUE.
  #' @param remove_asterisk Boolean. If TRUE, asterisks in 'ClustName' are removed. Default is TRUE.
  #' @examples cleanup_gencontext(prot, domains_rename, T, F)

  # Create cleaned up GenContext column
  prot$GenContext <- prot$GenContext.orig

  if (nrow(domains_rename) != 0) {
    ## Domains_rename
    for (x in 1:length(domains_rename$old)) {
      target <- domains_rename$old[x]
      replacement <- domains_rename$new[x]
      prot$GenContext <- prot$GenContext %>%
        str_replace_all(coll(target, TRUE), replacement)
    }
  }

  ## Reverse operons | Straighten them out!
  prot <- reverse_operon(prot)

  prot <- replaceQMs(prot, "GenContext")
  ## Optional parameters
  # Condense repeats
  if (repeat2s) {
    prot <- repeat2s(prot, "GenContext")
  }

  # Remove the Asterisks
  if (remove_asterisk) {
    prot <- remove_astrk(prot, colname = "GenContext")
  }

  return(prot)
}

cleanup_GeneDesc <- function(prot, column) {
  #' Return trailing period that occurs in GeneDesc column
  prot[, "GeneDesc"] <- gsub("\\.$", "", prot %>% pull(column))
  prot[, "GeneDesc"] <- gsub("%2C", ",", prot %>% pull(column))
  return(prot)
}


pick_longer_duplicate <- function(prot, column) {
  col <- sym(column)

  prot$row.orig <- 1:nrow(prot)

  # Get list of duplicates
  dups <- prot %>%
    group_by(AccNum) %>%
    summarize("count" = n()) %>%
    filter(count > 1) %>%
    arrange(-count) %>%
    merge(prot, by = "AccNum")

  dup_acc <- dups$AccNum

  longest_rows <- c()
  remove_rows <- c()
  for (acc in dup_acc) {
    dup_rows <- dups %>% filter(AccNum == acc)

    longest <- dup_rows[which(nchar(pull(dup_rows, {{ col }})) == max(nchar(pull(dup_rows, {{ col }}))))[1], "row.orig"]

    longest_rows <- c(longest_rows, longest)

    to_remove <- dup_rows[which(dup_rows$row.orig != longest), "row.orig"][]

    # dup_rows[which(nchar(pull(dup_rows,{{col}})) == max(nchar(pull(dup_rows,{{col}}))))[2:nrow(dup_rows)], "row.orig"]
    remove_rows <- c(remove_rows, to_remove)
  }

  # grab all the longest rows
  unique_dups <- prot[-remove_rows, ] %>% select(-row.orig)

  return(unique_dups)
}

cleanup_lineage <- function(prot, lins_rename) {
  for (i in 1:nrow(lins_rename)) {
    prot$Lineage <- gsub(lins_rename$old[i], lins_rename$new[i],
      x = prot$Lineage,
      ignore.case = T
    )
    # prot$Lineage = gsub(">REMOVE", "",
    #                           x = prot$Lineage,
    #                           ignore.case = T)
  }

  return(prot)
}
