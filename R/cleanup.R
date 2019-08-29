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

###########################
#### CLEANUP FUNCTIONS ####
###########################
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
cleanup_species <- function(prot, remove_empty = FALSE){
  # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  prot$Species <- prot$Species %>%
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
}

##########################
#'Cleanup Genomic Contexts
#'
#'Cleans up the GenContext column of a data frame by removing certain characters and rows.
#'
#'This function removes empty rows based on the 'GenContext' column.
#'A cleaned up version of the data table is returned.
#'
#'@param prot A data frame that contains columns 'GenContext'
#'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'GenContext' are removed
#'@examples cleanup_species(pspa.sub,TRUE)
cleanup_gencontext <- function(prot, remove_empty=TRUE){
  # prot$GenContext <- prot$GenContext.orig
  # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
  if(remove_empty){
    prot <- prot %>%
      as_tibble() %>%
      filter(grepl("\\*", GenContext)) %>%		# Keep only rows with Query (*)
      filter(!grepl("^-$", GenContext)) %>%		# remove "-"
      filter(!grepl("^NA$", GenContext)) %>%	# remove "NA"
      filter(!grepl("^$", GenContext)) #%>%		# remove empty rows
  }
return(prot)
}
#Switch case for remove.empty.rows, check efficiency
#Don't call other psp functions within these functions

###########################
#'Condense repeated domains
#'
#'Condenses repeated domains in the specified column.
#'
#'This function ...
#'Certain domains can be removed according to an additional data frame.
#'The original data frame is returned with the corresponding cleaned up column.
#'
#'@param prot A data frame containing a 'DomArch' column
#'@param by_column Can take values "DomArch" or "ClustName" to ...
#'@examples repeat2s(all, "DomArch")
repeat2s <- function(prot, by_column="DomArch"){
  switch(by_column,
         DomArch = (prot$DomArch <- prot$DomArch %>%
                      str_replace_all("\\+", " ") %>%
                      str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
                      str_replace_all(" ", "+")),
         ClustName = (prot$ClustName <- prot$ClustName %>%
                        str_replace_all("\\+", " ") %>%
                        str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
                        str_replace_all(" ", "+")))
  return(prot)
}

##############################
#'Cleanup Domain Architectures
#'
#'Cleans the DomArch column by replacing/removing certain domains
#'
#'This function cleans the DomArch column of one data frame by renaming certain domains according to a second data frame.
#'Certain domains can be removed according to an additional data frame.
#'The original data frame is returned with the clean DomArchs column and the old domains in the DomArchs.old column.
#'
#'@param prot A data frame containing a 'DomArch' column
#'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the corresponding replacement
#'values in a column 'new'
#'@param domains_remove A data frame containing the domain names to be removed in a column called 'domains'
#'@examples replace_doms(pspa.sub,domains.replace,domains.remove)
cleanup_domarch <- function(prot, domains_rename, domains_remove){ # was "replace_doms"
  prot$DomArch <- prot$DomArch.orig
  # Calling repeat2s
  repeat2s(prot=prot, by_column="DomArch")
  #replace domains based on the domains_rename list
  for(j in 1:length(domains_rename$old)){
    prot$DomArch <- str_replace_all(prot$DomArch,
                                    as.vector(domains_rename$old[j]),
                                    as.vector(domains_rename$new[j]))
  }
  #remove domains based on the domains_remove list
  for(j in 1:length(as.vector(domains_remove$domains))){
    prot$DomArch <- str_remove_all(prot$DomArch,
                                   as.vector(domains_remove$domains[j]))
  }
  #remove '+' at the start and end, as well as consecuative '+'
  prot$DomArch <- gsub("\\++\\+","\\+", prot$DomArch)
  prot$DomArch <- gsub("^\\+","", prot$DomArch)
  prot$DomArch <- gsub("\\+$","", prot$DomArch)

  return(prot)
  # return(cbind(prot, DomArch.old))
}

######################
#'Cleanup cluster file
#'
#'Cleans a cluster file by removing rows that do not contain the query in the cluster.
#'
#'This function removes irrelevant rows which do not contain the query protein within the ClustName column.
#'The return value is the cleaned up data frame.
#'
#'@param cls_data A data frame that must contain columns Query and ClustName.
#'@param by_column Can take values "ClustName" or "DomArch" to filter either cluster names
#' or domain architectures by the query name! Default: ClustName.
#'@examples cleanup_clust(all, "ClustName")
cleanup_clust <- function(cls_data, by_column="ClustName"){
  master <- cls_data[0,]
  #colnames(master) <- colnames(cls_data)

  # retrieve unique query names
  queryNames <- select(cls_data, Query) %>%
    distinct()

  for(i in 1:length(queryNames$Query)){
    #print(queryNames[i,])
    filtered_by_query <- cls_data %>%
      filter(Query == as.character(queryNames[i, "Query"]))

    switch(by_column,
           # Return BLASTCLUST names that match the query domain
           ClustName = (master <- rbind(master,
                                        filtered_by_query[grep(queryNames[i,],
                                                               filtered_by_query$ClustName,
                                                               ignore.case = T),])),
           # Return domain architectures that match the query domain
           DomArch = (master <- rbind(master,
                                      filtered_by_query[grep(queryNames[i,],
                                                             filtered_by_query$DomArch,
                                                             ignore.case = T),])))
  }

  # !!UNFIXED ISSUE!! SIG+TM+TM+... kind of architectures without explicit domain names are lost.
  # Need a way to take care of true hits that don't go by the same name.
  return(master)
}

