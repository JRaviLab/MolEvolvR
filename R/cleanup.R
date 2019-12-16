## Functions to clean up .op_ins_cls files by:
## Species, ClustName, DomArch, GenContext
## To create consistent names and take care of repeats & remove empty rows
## Created: Aug 11, 2017
## Modified: Dec 11, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")

###########################
#### CLEANUP FUNCTIONS ####
###########################
remove_empty <- function(prot, by_column="DomArch"){
  #'Remove empty rows by column
  #'
  #'Removes empty rows in the specified column.
  #'
  #'This function ...
  #'The original data frame is returned with the corresponding cleaned up column.
  #'
  #'@param prot A data frame containing 'DomArch', 'Species', 'GenContext', 'ClustName' columns.
  #'@param by_column Column by which empty rows should be removed to domain+domain -> domain(s).
  #' Default column is 'DomArch'. Can also take the following as input, 'Species', 'GenContext', 'ClustName'.
  #'@examples remove_empty(prot, "DomArch")

  #Switch case for remove_empty, check efficiency
  #Don't call other psp functions within these functions
  #! FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  #! Make a separate function of out of this?
  prot <- prot %>%
    as_tibble() %>%
    #filter(grepl("\\*", {{by_column}})) %>%		  # Keep only rows with Query (*) for GenContext
    filter(!grepl("^-$", {{by_column}})) %>%		# remove "-"
    filter(!grepl("^NA$", {{by_column}})) %>%	  # remove "NA"
    filter(!grepl("^$", {{by_column}}))       	# remove empty rows

  return(prot)
}

###########################
repeat2s <- function(prot, by_column="DomArch"){
  #'Condense repeated domains
  #'
  #'Condenses repeated domains in the specified column.
  #'
  #'This function ...
  #'Certain domains can be removed according to an additional data frame.
  #'The original data frame is returned with the corresponding cleaned up column.
  #'
  #'@param prot A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.
  #'@param by_column Column in which repeats are condensed to domain+domain -> domain(s).
  #' Default column is 'DomArch'. Can also take the following as input, 'GenContext', 'ClustName'.
  #'@examples repeat2s(prot, "DomArch")
  prot[,by_column] <- prot[,by_column] %>%
    mutate_all(funs(str_replace_all(.,
                                    pattern="\\+",
                                    replacement=" "))) %>%
    mutate_all(funs(str_replace_all(.,
                                    pattern="(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+",
                                    replacement="\\1(s)"))) %>%
    mutate_all(funs(str_replace_all(.,
                                    pattern=" ",
                                    replacement="+")))
  return(prot)
}

###########################
remove_tails <- function(prot, by_column="DomArch",
                         domains_keep=domains_keep){ #!! currently redundant
  #'Remove tails/singletons
  #'
  #'
  #'This function ...
  #'Certain low frequency domain architectures can be removed.
  #'The original data frame is returned with the corresponding cleaned up column.
  #'
  #'@param prot A data frame containing 'DomArch', 'GenContext', 'ClustName' columns.
  #'@param by_column Default column is 'DomArch'. Can also take 'ClustName', 'GenContext' as input.
  #'@param keep_query Default is TRUE. Keeps tail entries that contain the query domains.
  #'@examples remove_tails(prot, "DomArch")
  by_column <- sym(by_column)
  domain_count <- prot %>%
    group_by({{by_column}}) %>%
    summarize(count = n())

  ## Identify tails
  tails <- domain_count %>% filter(count == 1)

  ## Domains_keep
  # Keep tails with query domains
  # Contains all domains separated by "|"
  domains_for_grep <- paste(domains_keep$domains, collapse = "|")
  # Remove rows with no domains contained within domains_keep
  tails <- tails %>%
    filter(!grepl(domains_for_grep, {{by_column}})) ## CRAZY thing doesn't work!

  # Remove tails
  prot <- prot %>% filter(!({{by_column}} %in% tails))

  return(prot)
}

###########################
cleanup_species <- function(prot, remove_empty=FALSE){
  #'Cleanup Species
  #'
  #'Cleans up the species column of a data frame by removing certain characters and rows.
  #'
  #'This function removes unneccessary characters from the 'Species' column.
  #'A cleaned up version of the data table is returned.
  #'
  #'@param prot A data frame that contains columns 'Species'.
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'Species' are removed.
  #' Default is false.
  #'@examples cleanup_species(prot,TRUE)
  # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  prot$Species <- prot$Species.orig %>%
    str_replace_all(coll("sp. ", TRUE), "sp ") %>%
    str_replace_all(coll("str. ", TRUE), "str ") %>%
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
    str_replace_all(coll("  ", TRUE), " ")

  # ADD SOMETHING (!) for removing empty rows
  # !! CHECK !! Species vs Species_old
  if(remove_empty){
    prot <- remove_empty(prot=prot, by_column="Species")
  }

  return(prot)
}

######################
cleanup_clust <- function(prot,
                          domains_keep, domains_rename,
                          repeat2s=TRUE, remove_tails = FALSE,
                          remove_empty=FALSE){
  #'Cleanup cluster file
  #'
  #'Cleans a cluster file by removing rows that do not contain the query in the cluster.
  #'
  #'This function removes irrelevant rows which do not contain the query protein within the ClustName column.
  #'The return value is the cleaned up data frame.
  #'
  #'@param prot A data frame that must contain columns Query and ClustName.
  #'@param domains_keep A data frame containing the domain names to be retained.
  #'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the
  #'corresponding replacement values in a column 'new'.
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'ClustName' are condensed. Default is TRUE.
  #'@param remove_tails Boolean. If TRUE, 'ClustName' will be filtered based on domains to keep/remove. Default is FALSE.
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'ClustName' are removed. Default is FALSE.
  #'@examples cleanup_clust(prot, TRUE, FALSE, domains_keep, domains_rename)


  prot$ClustName <- prot$ClustName.orig

  ## Basic Cleanup
  # Remove '+' at the start and end, as well as consecuative '+'
  prot$ClustName <- gsub("\\++\\+","\\+", prot$ClustName)
  prot$ClustName <- gsub("^\\+","", prot$ClustName)
  prot$ClustName <- gsub("\\+$","", prot$ClustName)

  ## Domains_rename
  for(x in 1:length(domains_rename$old)){
    target <- domains_rename$old[x]
    replacement <- domains_rename$new[x]
    prot$ClustName <- prot$ClustName %>% str_replace_all(target,replacement)
  }

  ## Domains_keep
  # Character for greping for rows with domains_keep
  # Contains all domains separated by "|"
  domains_for_grep <- paste(domains_keep$domains, collapse = "|")
  # Remove rows with no domains contained within domains_keep
  prot <- prot %>% filter(grepl(domains_for_grep,ClustName))

  ## Optional parameters
  # Condense repeats
  if(repeat2s){
    prot <- repeat2s(prot, by_column="ClustName")
  }
  # Remove singletons
  #if(remove_tails){
  #  prot <- prot %>% filter(!grepl(".1$", ClustID))
  #}
  if(remove_tails){
    prot <- remove_tails(prot, by_column="ClustName")
  }
  # Remove empty rows
  if(remove_empty){
    prot <- remove_empty(prot=prot, by_column="ClustName")
  }

  # !!UNFIXED ISSUE!! Currently requires manual intervention!
  # SIG+TM+TM+... kind of architectures without explicit domain names are lost.
  # Need a way to take care of true hits that don't go by the expected domain name.
  return(prot)
}

##############################
cleanup_domarch <- function(prot,
                            domains_keep, domains_rename,
                            repeat2s=TRUE, remove_tails = FALSE,
                            remove_empty=F,
                            domains_ignore=NULL){
  #'Cleanup Domain Architectures
  #'
  #'Cleans the DomArch column by replacing/removing certain domains
  #'
  #'This function cleans the DomArch column of one data frame by renaming certain domains according to a second data frame.
  #'Certain domains can be removed according to an additional data frame.
  #'The original data frame is returned with the clean DomArchs column and the old domains in the DomArchs.old column.
  #'
  #'@param prot A data frame containing a 'DomArch' column
  #'@param domains_keep A data frame containing the domain names to be retained.
  #'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the
  #'corresponding replacement values in a column 'new'.
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'DomArch' are condensed. Default is TRUE.
  #'@param remove_tails Boolean. If TRUE, 'ClustName' will be filtered based on domains to keep/remove. Default is FALSE.
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'DomArch' are removed. Default is FALSE.
  #'@param domains_ignore A data frame containing the domain names to be removed in a column called 'domains'
  #'@examples cleanup_domarch(prot, TRUE, FALSE, domains_keep, domains_rename, domains_ignore=NULL)

  prot$DomArch <- prot$DomArch.orig

  ## Basic Cleanup
  # Remove '+' at the start and end, as well as consecuative '+'
  prot$DomArch <- gsub("\\++\\+","\\+", prot$DomArch)
  prot$DomArch <- gsub("^\\+","", prot$DomArch)
  prot$DomArch <- gsub("\\+$","", prot$DomArch)

  ## Domains_rename
  # Replace domains based on the domains_rename list
  for(j in 1:length(domains_rename$old)){
    prot$DomArch <- str_replace_all(prot$DomArch,
                                    as.vector(domains_rename$old[j]),
                                    as.vector(domains_rename$new[j]))
  }

  ## Domains_keep
  # Character for greping for rows with domains_keep
  # Contains all domains separated by "|"
  domains_for_grep <- paste(domains_keep$domains, collapse = "|")
  # Remove rows with no domains contained within domains_keep
  prot <- prot %>% filter(grepl(domains_for_grep, DomArch))

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
  if(remove_tails){
    prot <- remove_tails(prot=prot, by_column="DomArch")
  }
  # Condense repeats
  if(repeat2s){
    ## Error in UseMethod("tbl_vars") : no applicable method for 'tbl_vars' applied to an object of class "character"
    prot <- repeat2s(prot=prot, by_column="DomArch")
  }
  # Remove empty rows
  #! FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  #! Make a separate function of out of this?
  if(remove_empty){
    prot <- remove_empty(prot=prot, by_column="DomArch")
  }

  return(prot)
}

##########################
cleanup_gencontext <- function(prot, domains_rename,
                               repeat2s=TRUE){
  #'Cleanup Genomic Contexts
  #'
  #'Cleans up the GenContext column of a data frame by removing certain characters and rows.
  #'
  #'This function removes empty rows based on the 'GenContext' column.
  #'A cleaned up version of the data table is returned.
  #'
  #'@param prot A data frame that contains columns 'GenContext.orig'
  #'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'GenContext' are condensed. Default is TRUE.
  #'@examples cleanup_gencontext(prot, domains_rename, T, F)
  prot$GenContext <- prot$GenContext.orig

  ## Domains_rename
  for(x in 1:length(domains_rename$old)){
    target <- domains_rename$old[x]
    replacement <- domains_rename$new[x]
    prot$GenContext <- prot$GenContext %>% str_replace_all(target,replacement)
  }

  ## Reverse operons | Straighten them out!
  prot <- reverse_operon(prot)

  ## Optional parameters
  # Condense repeats
  if(repeat2s){
    prot <- repeat2s(prot, "GenContext")
  }

  return(prot)
}
