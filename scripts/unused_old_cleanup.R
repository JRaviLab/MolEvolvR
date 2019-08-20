# old names for "op_ins_cls" files
colnames.op_ins_cls <- function(){
  c("AccNum", "GenContext",
    "SIG.TM.PFAM", "SIG.TM.LADB", "PFAM",
    "Length", "GenName", "Lineage", "Species",
    "Annotation", "GI")
}

## Function to replace singleton domains with the (s) form
domarch.convert2s.forupset <- function(x){ x %>%
    str_replace_all("^TM\\+","TM(s)+") %>%
    str_replace_all("\\+TM\\+","+TM(s)+") %>%
    str_replace_all("TM$","TM(s)") # %>%
  # str_replace_all("PspC\\+","PspC(s)+") %>%
  # str_replace_all("PspC$","PspC(s)")
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


##################################

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