## Functions to clean up .op_ins_cls files
## To create consistent names and take care of repeats
## To get element and word counts

## Modified: Jun 06, 2019 (some func have been moved to summarize.R)
## Created: Aug 11, 2017
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")

###########################
#### CLEANUP FUNCTIONS ####
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
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'Species' are removed
  #'@examples cleanup_species(pspa,TRUE)
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

  # ADD SOMETHING for removing empty rows
  # !! CHECK !! Species vs Species_old
  if(remove_empty){
    prot <- prot %>%
      as_tibble() %>%
      filter(!grepl("^-$", Species)) %>%		# remove "-"
      filter(!grepl("^NA$", Species)) %>%	# remove "NA"
      filter(!grepl("^$", Species)) #%>%		# remove empty rows
  }
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

##########################
cleanup_gencontext <- function(prot, repeat2s=TRUE, remove_empty=FALSE, domains_rename){
  #'Cleanup Genomic Contexts
  #'
  #'Cleans up the GenContext column of a data frame by removing certain characters and rows.
  #'
  #'This function removes empty rows based on the 'GenContext' column.
  #'A cleaned up version of the data table is returned.
  #'
  #'@param prot A data frame that contains columns 'GenContext.orig'
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'GenContext' are condensed. Default is TRUE.
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'GenContext' are removed
  #'@examples cleanup_gencontext(prot, FALSE)
  prot$GenContext <- prot$GenContext.orig
  # Condsense repeatsd
  if(repeat2s){
    prot <- repeat2s(prot, "GenContext")
  }

  for(x in 1:length(domains_rename$old)){
    target <- domains_rename$old[x]
    replacement <- domains_rename$new[x]
    prot$GenContext <- prot$GenContext %>% str_replace_all(target,replacement)
  }

  # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  if(remove_empty){
    prot <- prot %>%
      as_tibble() %>%
      filter(grepl("\\*", GenContext)) %>%		# Keep only rows with Query (*)
      filter(!grepl("^-$", GenContext)) %>%		# remove "-"
      filter(!grepl("^NA$", GenContext)) %>%	# remove "NA"
      filter(!grepl("^$", GenContext)) #%>%		# remove empty rows
  }

  prot <- reverse_operon(prot)
  return(prot)
}
#Switch case for remove.empty.rows, check efficiency
#Don't call other psp functions within these functions

##############################
cleanup_domarch <- function(prot, repeat2s=TRUE, domains_rename, domains_ignore = NULL){
  #'Cleanup Domain Architectures
  #'
  #'Cleans the DomArch column by replacing/removing certain domains
  #'
  #'This function cleans the DomArch column of one data frame by renaming certain domains according to a second data frame.
  #'Certain domains can be removed according to an additional data frame.
  #'The original data frame is returned with the clean DomArchs column and the old domains in the DomArchs.old column.
  #'
  #'@param prot A data frame containing a 'DomArch' column
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'DomArch' are condensed. Default is TRUE.
  #'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the corresponding replacement
  #'values in a column 'new'
  #'@param domains_ignore A data frame containing the domain names to be removed in a column called 'domains'
  #'@examples cleanup_domarch(prot, TRUE, domains_rename, domains_ignore)

  prot$DomArch <- prot$DomArch.orig

  # Condense repeats
  if(repeat2s){
    ## Error in UseMethod("tbl_vars") : no applicable method for 'tbl_vars' applied to an object of class "character"
    prot <- repeat2s(prot=prot, by_column="DomArch")
  }

  # Replace domains based on the domains_rename list
  for(j in 1:length(domains_rename$old)){
    prot$DomArch <- str_replace_all(prot$DomArch,
                                    as.vector(domains_rename$old[j]),
                                    as.vector(domains_rename$new[j]))
  }

  # Remove domains based on the domains_ignore list
  # ?? Remove the domains or the rows? Check, please!
  if( !is.null(domains_ignore)){
    for(j in 1:length(domains_ignore$domains)){
      prot$DomArch <- str_remove_all(prot$DomArch,
                                     as.vector(domains_ignore$domains[j]))
    }
  }
  # Remove '+' at the start and end, as well as consecuative '+'
  prot$DomArch <- gsub("\\++\\+","\\+", prot$DomArch)
  prot$DomArch <- gsub("^\\+","", prot$DomArch)
  prot$DomArch <- gsub("\\+$","", prot$DomArch)

  return(prot)
}

######################
cleanup_clust <- function(cls_data, repeat2s=TRUE,
                          #remove_tails = FALSE,
                          domains_keep, domains_rename){
  #'Cleanup cluster file
  #'
  #'Cleans a cluster file by removing rows that do not contain the query in the cluster.
  #'
  #'This function removes irrelevant rows which do not contain the query protein within the ClustName column.
  #'The return value is the cleaned up data frame.
  #'
  #'@param cls_data A data frame that must contain columns Query and ClustName.
  #'@param by_column Boolean. If TRUE, 'ClustName' will be filtered based on domains to keep/remove. Default is TRUE.
  #'@param repeat2s Boolean. If TRUE, repeated domains in 'ClustName' are condensed. Default is TRUE.
  #'@examples cleanup_clust(prot, "ClustName", TRUE)


  cls_data$ClustName <- cls_data$ClustName.orig

  # Character for greping for rows with domains_keep
  # Contains all domains separated by "|"
  domains_for_grep <- paste(domains_keep$domains, collapse = "|")

  # Remove rows with no domains contained within domains_keep
  cls_data <- cls_data %>% filter(grepl(domains_for_grep,ClustName))


  for(x in 1:length(domains_rename$old)){
    target <- domains_rename$old[x]
    replacement <- domains_rename$new[x]
    cls_data$ClustName <- cls_data$ClustName %>% str_replace_all(target,replacement)
  }


  # Condense repeats
  if(repeat2s){
    cls_data <- repeat2s(cls_data, "ClustName")
  }

  #if(remove_tails){
  #  cls_data <- cls_data %>% filter(!grepl(".1$", ClustID))
  #}

  # !!UNFIXED ISSUE!! SIG+TM+TM+... kind of architectures without explicit domain names are lost.
  # Need a way to take care of true hits that don't go by the expected domain name.
  return(cls_data)
}


remove_tails <- function(prot){
  domain_count <- prot %>% group_by(DomArch) %>% summarize(count = n())
  tails <- domain_count %>% filter(count == 1)
  prot <- prot %>% filter(!(DomArch %in% tails))
}