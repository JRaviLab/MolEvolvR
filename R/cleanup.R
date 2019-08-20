## Functions to clean up .op_ins_cls files
## To create consistent names and take care of repeats
## To get element and word counts

## Modified: Jun 06, 2019 (some func have been moved to summarize.R)
## Created: Aug 11, 2017
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzchen)

###########################
#### CLEANUP FUNCTIONS ####
###########################
#'Cleanup Species
#'
#'Cleans up the species column of a data frame by removing certain characters and rows.
#'
#'This function removes unneccessary characters from the 'Species' column.
#'Certain rows may also removed from the table based on values in the 'GenContext.norep' column.
#'A cleaned up version of the data table is returned.
#'
#'@param prot A data frame that contains columns 'Species' and 'GenContext.norep'
#'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'GenContext.norep' are removed
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

  ##########################
  #'Cleanup Genomic Contexts
  #'
  #'Cleans up the GenContext column of a data frame by removing certain characters and rows.
  #'
  #'This function removes empty rows based on the 'GenContext' column.
  #'A cleaned up version of the data table is returned.
  #'
  #'@param prot A data frame that contains columns 'GenContext.norep'
  #'@param remove_empty Boolean. If TRUE, rows with empty/unnecessary values in 'GenContext.norep' are removed
  #'@examples cleanup_species(pspa.sub,TRUE)
  cleanup_gencontext <- function(prot, remove_empty=TRUE)
    # FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
    if(remove_empty){
      prot <- prot %>%
        as_tibble() %>%
        filter(grepl("\\*", GenContext.norep)) %>%		# Keep only rows with Query (*)
        filter(!grepl("^-$", GenContext.norep)) %>%		# remove "-"
        filter(!grepl("^NA$", GenContext.norep)) %>%	# remove "NA"
        filter(!grepl("^$", GenContext.norep)) #%>%		# remove empty rows
    }
  return(prot)
}
#Switch case for remove.empty.rows, check efficiency
#Don't call other psp functions within these functions

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
cleanup_domarch <- function(prot,domains_rename, domains_remove){ # was "replace_doms"
  DomArch.old <- prot$DomArch

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

  # Deal with REPEATED DOMAINS by replacing with (s)!!
  prot$DomArch <- prot$DomArch %>%
    str_replace_all("\\+", " ") %>%
    str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
    str_replace_all(" ", "+")


  return(cbind(prot, DomArch.old))
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
#'
#'@examples cls_cleanup(all_op_ins)
cleanup_clust <- function(cls_data){
  master <- cls_data[0,]
  #colnames(master) = colnames(cls_data)
  queryNames <- select(cls_data, Query) %>%
    distinct()

  for(i in 1:length(queryNames$Query)){
    #print(queryNames[i,])
    filtered_by_query <- cls_data %>%
      filter(Query == as.character(queryNames[i, "Query"]))
    master <- rbind(master,
                    filtered_by_query[grep(queryNames[i,],
                                           filtered_by_query$ClustName,
                                           ignore.case = T),])
  }

  return(master)
}

