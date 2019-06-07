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
species.cleanup <- function(x) {
	x %>%
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
}
## Function to cleanup Toast-rack names
# replace.toastrack <- function(x) {
# 	x %>%
# 		str_replace_all(coll("toast-rack", TRUE), "Toast-rack") %>%
# 		str_replace_all(coll("toastrack", TRUE), "Toast-rack") %>%
# 		str_replace_all("DUF4097", "Toast-rack") %>%
# 		str_replace_all("DUF2154", "Toast-rack") %>%
# 		str_replace_all("DUF2807", "Toast-rack") %>%
# 		str_replace_all("DUF2157", "Toast_rack_N")
# }

## Function to deal with REPEATED DOMAINS with (s)!!!
repeats2s <- function(x){
	x %>%
		str_replace_all("\\+", " ") %>%
		str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
		str_replace_all(" ", "+")
}
## Function to remove spurious domains (from ignore.list) | DomArch column !!!
## ALERT !!! See if this creates problems!!!
remove.toastrack.ignored.doms <- function(x) {
	x %>%
		str_replace_all("EP1\\+","") %>%
		str_replace_all("EIID-AGA\\+","") %>%
		str_replace_all("LPD29\\+","") %>%
		str_replace_all("Imm3\\+","") %>%
		# str_replace_all("FilH","") %>%
		str_replace_all("ASH-IG", "IG")
}
## Function to replace singleton domains with the (s) form
domarch.convert2s.forupset <- function(x){ x %>%
		str_replace_all("^TM\\+","TM(s)+") %>%
		str_replace_all("\\+TM\\+","+TM(s)+") %>%
		str_replace_all("TM$","TM(s)") # %>%
		# str_replace_all("PspC\\+","PspC(s)+") %>%
		# str_replace_all("PspC$","PspC(s)")
}
## Function to remove empty rows
remove.empty.rows <- function(x){	x %>%
		as_tibble() %>%
		filter(grepl("\\*", GenContext.norep)) %>%		# Keep only rows with Query*
		filter(!grepl("^-$", GenContext.norep)) %>%		# remove "-"
		filter(!grepl("^NA$", GenContext.norep)) %>%	# remove "NA"
		filter(!grepl("^$", GenContext.norep)) #%>%		# remove empty rows
		# filter(!grepl("^-$", DomArch.norep)) %>%			# remove "-"
		# filter(!grepl("^NA$", DomArch.norep)) %>%			# remove "NA"
		# filter(!grepl("^$", DomArch.norep))						# remove empty rows
}




#Function that applies both remove.empty.rows and species.cleanup to a data table
cleanup_species <- function(prot){
  prot$Species <- species.cleanup(prot$Species)
  return(remove.empty.rows(prot))

  # #species.cleanup
  # prot$Species %>%
  #   str_replace_all(coll("sp. ", TRUE), "sp ") %>%
  #   str_replace_all(coll("str. ", TRUE), "str ") %>%
  #   str_replace_all(coll(" = ", TRUE), " ") %>%
  #   str_replace_all(coll("-", TRUE), "") %>%
  #   str_replace_all(coll(".", TRUE), "") %>%
  #   str_replace_all(coll("(", TRUE), "") %>%
  #   str_replace_all(coll(")", TRUE), "") %>%
  #   str_replace_all(coll("[", TRUE), "") %>%
  #   str_replace_all(coll("]", TRUE), "") %>%
  #   str_replace_all(coll("\\", TRUE), "") %>%
  #   str_replace_all(coll("/", TRUE), "") %>%
  #   str_replace_all(coll("\'", TRUE), "") %>%
  #   str_replace_all(coll("  ", TRUE), " ")
  #   #Remove.empty.rows
  # prot %>%
  #   as_tibble() %>%
  #   filter(grepl("\\*", GenContext.norep)) %>%		# Keep only rows with Query*
  #   filter(!grepl("^-$", GenContext.norep)) %>%		# remove "-"
  #   filter(!grepl("^NA$", GenContext.norep)) %>%	# remove "NA"
  #   filter(!grepl("^$", GenContext.norep)) #%>%		# remove empty rows
  # return(prot)
}


#'Cleaning Domain Architecture
#'
#'Cleans the DomArch column by replacing/removing certain domains
#'
#'This function cleans the DomArch column of one data frame by renaming certain domains according to a second data frame.
#'Additionally, repeated domains are condensed into the name of the domain with (s) appended. The original data
#'frame is returned with the clean DomArchs column and the old domains in the DomArchs.old column.
#'
#'@param eDNA_data A data frame containing a 'DomArch' column
#'@param domains_rename A data frame containing the domain names to be replaced in a column 'old' and the corresponding replacement
#'values in a column 'new'
#'@examples replace_doms(pspa.sub,pspdomains.rename)
replace_doms <- function(eDNA_data,domains_rename){
  DomArch.old <- eDNA_data$DomArch
  for (j in 1:length(domains_rename$old)) {

    #Find and replace the occurances of domain where it starts with domain+
    regterm=paste ("^",domains_rename$old[j],"\\+",sep="")
    #repterm is the string "domofnew+"
    repterm=paste (domains_rename$new[j],"\\+",sep="")
    #All indeces in eDNA_data$arch that start with regterm, are replaced with repterm
    eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)]=gsub(regterm, repterm,x = eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)])

    #Find and replace the occurances of domain where it ends in +domain
    regterm=paste ("\\+",domains_rename$old[j],"$",sep="")
    #repterm is now "+domainofnew" (very similar to first repterm)
    repterm=paste ("\\+",domains_rename$new[j],sep="")
    #All indeces in eDNA_data$arch that match second regterm are replaced with second repterm
    eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)]=gsub(regterm, repterm,x = eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)])

    #What if domain is replaced by "" -- there would be nothing there
    #regterm3 = "+domainofold+"
    regterm=paste ("\\+",domains_rename$old[j],"\\+",sep="")
    repterm=paste ("\\+",domains_rename$new[j],"\\+",sep="")
    eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)]=gsub(regterm, repterm,x = eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)])

    # #all occurances where domain to be renamed occurs by itself is replaced
    # regterm=paste (domains_rename$old[j],sep="")
    # repterm=paste (domains_rename$new[j],sep="")
    # eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)]=gsub(regterm, repterm,x = eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)])

    #repterm= domains_rename$new[j]
    #What does[1] do?
    #While loop b/c it can appear multiple times
    #while iterates keeps going until no more strings in eDNA_data$arch match regterm
    while ( length(eDNA_data$DomArch[grep(regterm[1], eDNA_data$DomArch)]) > 0) {
      #replace the regterm with repterm
      eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)]=gsub(regterm, repterm,x = eDNA_data$DomArch[grep(regterm, eDNA_data$DomArch)])

    }

  }

  eDNA_data$DomArch <- repeats2s(eDNA_data$DomArch) %>%
    domarch.convert2s.forupset()
  #eDNA_data$DomArch <- domarch.convert2s.forupset(eDNA_data$DomArch)
  return(cbind(eDNA_data,DomArch.old))
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


