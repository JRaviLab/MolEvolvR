## Modified: Apr 08, 2019
## Created: Aug 11, 2017
## Janani Ravi (@jananiravi)
## Functions to clean up .op_ins_cls files
## To create consistent names and take care of repeats
## To get element and word counts

##############
## COLNAMES ##
##############
## FUNCTION to ASSIGN COLUMN NAMES based on AUG 2017 VA format
colnames.op_ins_cls <- function(){
	c("AccNum", "GenContext",
		"SIG.TM.PFAM", "SIG.TM.LADB", "PFAM",
		"Length", "GenName", "Lineage", "Species",
		"Annotation", "GI")
}

###########################
#### CLEANUP FUNCTIONS ####
###########################
## Function to cleanup Species
# species.cleanup <- function(x) {
# 	x %>%
# 		str_replace_all(coll("sp. ", TRUE), "sp ") %>%
# 		str_replace_all(coll("str. ", TRUE), "str ") %>%
# 		str_replace_all(coll(" = ", TRUE), " ") %>%
# 		str_replace_all(coll("-", TRUE), "") %>%
# 		str_replace_all(coll(".", TRUE), "") %>%
# 		str_replace_all(coll("(", TRUE), "") %>%
# 		str_replace_all(coll(")", TRUE), "") %>%
# 		str_replace_all(coll("[", TRUE), "") %>%
# 		str_replace_all(coll("]", TRUE), "") %>%
# 		str_replace_all(coll("\\", TRUE), "") %>%
# 		str_replace_all(coll("/", TRUE), "") %>%
# 		str_replace_all(coll("\'", TRUE), "") %>%
# 		str_replace_all(coll("  ", TRUE), " ")
# }
## Function to cleanup Toast-rack names
##
# replace.toastrack <- function(x) {
# 	x %>%
# 		str_replace_all(coll("toast-rack", TRUE), "Toast-rack") %>%
# 		str_replace_all(coll("toastrack", TRUE), "Toast-rack") %>%
# 		str_replace_all("DUF4097", "Toast-rack") %>%
# 		str_replace_all("DUF2154", "Toast-rack") %>%
# 		str_replace_all("DUF2807", "Toast-rack") %>%
# 		str_replace_all("DUF2157", "Toast_rack_N")
# str_replace_all("ASH-IG", "IG")
# }

#use in
## Function to deal with REPEATED DOMAINS with (s)!!!
# repeats2s <- function(x){
# 	x %>%
# 		str_replace_all("\\+", " ") %>%
# 		str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
# 		str_replace_all(" ", "+")
# }

## Remove extra pluses that may occur from removal
##New file
## Function to remove spurious domains (from ignore.list) | DomArch column !!!
## ALERT !!! See if this creates problems!!!
# remove.toastrack.ignored.doms <- function(x) {
# 	x %>%
# 		str_replace_all("EP1\\+","") %>%
# 		str_replace_all("EIID-AGA\\+","") %>%
# 		str_replace_all("LPD29\\+","") %>%
# 		str_replace_all("Imm3\\+","")
#     #%>%
# 		# str_replace_all("FilH","") %>%
#
# }


## Function to replace singleton domains with the (s) form
domarch.convert2s.forupset <- function(x){ x %>%
		str_replace_all("^TM\\+","TM(s)+") %>%
		str_replace_all("\\+TM\\+","+TM(s)+") %>%
		str_replace_all("TM$","TM(s)") # %>%
		# str_replace_all("PspC\\+","PspC(s)+") %>%
		# str_replace_all("PspC$","PspC(s)")
}
## Function to remove empty rows
# remove.empty.rows <- function(x){	x %>%
# 		as_tibble() %>%
# 		filter(grepl("\\*", GenContext.norep)) %>%		# Keep only rows with Query*
# 		filter(!grepl("^-$", GenContext.norep)) %>%		# remove "-"
# 		filter(!grepl("^NA$", GenContext.norep)) %>%	# remove "NA"
# 		filter(!grepl("^$", GenContext.norep)) #%>%		# remove empty rows
# 		# filter(!grepl("^-$", DomArch.norep)) %>%			# remove "-"
# 		# filter(!grepl("^NA$", DomArch.norep)) %>%			# remove "NA"
# 		# filter(!grepl("^$", DomArch.norep))						# remove empty rows
# }



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
#'@examples cleanup_species(pspa.sub,TRUE)
cleanup_species <- function(prot,remove_empty = FALSE){
  #FUNCTIONS CALLED HERE, if else might be better since only two options, T and F
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
  if(remove_empty){
    prot <- prot%>%
      as_tibble() %>%
      filter(grepl("\\*", GenContext.norep)) %>%		# Keep only rows with Query*
      filter(!grepl("^-$", GenContext.norep)) %>%		# remove "-"
      filter(!grepl("^NA$", GenContext.norep)) %>%	# remove "NA"
      filter(!grepl("^$", GenContext.norep)) #%>%		# remove empty rows
  }
  return(prot)
}
#Switch case for remove.empty.rows, check efficiency
#Don't call other psp functions within these functions

#'Cleaning Domain Architecture
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
replace_doms <- function(prot,domains_rename, domains_remove){
  DomArch.old <- prot$DomArch

  #replace domains based on the domains_rename list
  for(j in 1:length(domains_rename$old)){
    prot$DomArch <- str_replace_all(prot$DomArch,as.vector(domains_rename$old[j]),as.vector(domains_rename$new[j]))
  }
  #remove domains based on the domains_remove list
  for(j in 1:length(as.vector(domains_remove$domains))){
    prot$DomArch <- str_remove_all(prot$DomArch,as.vector(domains_remove$domains[j]))
  }
  #remove '+' at the start and end, as well as consecuative '+'
  prot$DomArch <- gsub("\\++\\+","\\+",prot$DomArch)
  prot$DomArch <- gsub("^\\+","",prot$DomArch)
  prot$DomArch <- gsub("\\+$","",prot$DomArch)

  # Deal with REPEATED DOMAINS by replacing with (s)!!
  prot$DomArch <- prot$DomArch %>%
    str_replace_all("\\+", " ") %>%
    str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
    str_replace_all(" ", "+")


  return(cbind(prot,DomArch.old))
}


##################################
## Descriptions for functions ####
##################################
## colnames.op_ins_cls: FUNCTION to ASSIGN COLUMN NAMES based on AUG 2017 VA format
# cat("Colnames:
# AccNum, GenContext, SIG.TM.PFAM, SIG.TM.LADB, PFAM,
# Length, GenName, Lineage, Species, Annotation, GI")

#	## replace.toastrack: Function definition and calling | Prints for User
# 	cat("Renaming DUFs/TRs to Toast-rack & Toast_rack_N.
# \nFor e.g.:
# toast_rack$DomArch.norep %>%
# replace.toastrack()")

# 	## repeat2s: Function definition and calling | Prints for User
# 	cat("Converts repeats to (s).
# \nFor e.g.:
# query$SIG.TM.LADB %>%
#   repeat2s()
# \nquery$GenContext %>%
#   replace.toastrack() %>%  ## for toast_rack query
#   repeat2s()")

# ## remove.empty.rows: Function definition and calling | Prints for User
# cat("Removes empty rows from DomArch.norep & GenContext.norep columns.
# \nFor e.g.:
# query.sub <- query %>%
#   remove.empty.rows()")


