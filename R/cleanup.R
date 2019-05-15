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
replace.toastrack <- function(x) {
	x %>%
		str_replace_all(coll("toast-rack", TRUE), "Toast-rack") %>%
		str_replace_all(coll("toastrack", TRUE), "Toast-rack") %>%
		str_replace_all("DUF4097", "Toast-rack") %>%
		str_replace_all("DUF2154", "Toast-rack") %>%
		str_replace_all("DUF2807", "Toast-rack") %>%
		str_replace_all("DUF2157", "Toast_rack_N")
}

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
###########################
## COUNTS of DAs and GCs ##
## Before/after break up ##
###########################
## Function to obtain element counts (DA, GC)
counts_elements <- function(x, min.freq) {	x %>%
		table() %>%
		as_tibble() %>%
		`colnames<-`(c("elements", "freq")) %>%
		filter(!grepl("^-$", elements)) %>%		# remove "-"
		arrange(-freq) %>% filter(freq>=min.freq)
}

## Function to break up ELEMENTS to WORDS for DA and GC
elements2words <- function(x, type) {
	y <- x %>%
		str_replace_all("\\,"," ") %>%
		str_replace_all("\""," ")
	switch(type,
				 da2doms = {z <- y %>%
				 	str_replace_all("\\+"," ")},
				 gc2da = {z <- y %>%
				 	str_replace_all("\\<-"," ") %>%
				 	str_replace_all("-\\>"," ") %>%
				 	str_replace_all("\\|"," ")})
	# str_replace_all("^c\\($", " ") %>%		# remove "c("
	# str_replace_all("\\)$", " ") %>%			# remove ")"
	# str_replace_all("\\(s\\)"," ") %>%		# Ignoring repeats
	# str_replace_all("-"," ") %>%
	## replace \n, \r, \t
	z %>%
		str_replace_all("\n"," ") %>%
		str_replace_all("\r"," ") %>%
		str_replace_all("\t"," ") %>%
		## replace multiple spaces ...
		str_replace_all("    "," ") %>%
		str_replace_all("   "," ") %>%
		str_replace_all("  "," ") %>%
		str_replace_all("  "," ")
}

## Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
## to be used after elements2words
words2wc <- function(x){ x %>%
		str_replace_all("   "," ") %>%
		str_replace_all("  "," ") %>% str_replace_all("  "," ") %>%
		paste(collapse=" ") %>%
		strsplit(" ") %>%
		# filter(grepl(query.list[j], Query)) %>% # Create separate WCs for each Query
		# select(DA.wc) %>%
		table() %>% as_tibble() %>%
		`colnames<-`(c("words", "freq")) %>%
		## filter out 'spurious-looking' domains
		filter(!grepl(" \\{n\\}", words)) %>%
		filter(!grepl("^c\\($", words)) %>%		# remove "c("
		filter(!grepl("^\\)$", words)) %>%		# remove ")"
		filter(!grepl("^-$", words)) %>%			# remove "-"
		filter(!grepl("^$", words)) %>%				# remove empty rows
		filter(!grepl("^\\?$", words)) %>%		# remove "?"
		filter(!grepl("^\\?\\*$", words)) %>%	# remove "?*"
		filter(!grepl("^tRNA$", words)) %>%		# remove "tRNA"
		filter(!grepl("^ncRNA$", words)) %>%	# remove "ncRNA"
		filter(!grepl("^rRNA$", words)) %>%		# remove "rRNA"
		# filter(!grepl("\\*", words)) %>%			# Remove/Keep only Query
		arrange(-freq)
}
## Function to filter based on frequencies
filter.freq <- function(x, min.freq){ x %>%
		filter(freq>=min.freq)
}

#########################
## SUMMARY FUNCTIONS ####
## Changed Lineage to Lineage.final on Aug 31
#########################
## Function to summarize and retrieve counts by Domains & Domains+Lineage
summ.DA.byLin <- function(x) { x %>%
		filter(!grepl("^-$", DomArch.norep)) %>%
		group_by(DomArch.norep, Lineage.final) %>%
		summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
		arrange(desc(count))
}
summ.DA <- function(x){ x %>%
		group_by(DomArch.norep) %>%
		summarise(totalcount=sum(count), totallin=n()) %>% # totallin=n_distinct(Lineage),
		arrange(desc(totallin), desc(totalcount)) %>%
		filter(!grepl(" \\{n\\}",DomArch.norep)) %>%
		filter(!grepl("^-$", DomArch.norep))
}
summ.GC.byDALin <- function(x) { x %>%
		filter(!grepl("^-$", GenContext.norep)) %>%
		filter(!grepl("^-$", DomArch.norep)) %>%
		filter(!grepl("^-$", Lineage.final)) %>% filter(!grepl("^NA$", DomArch.norep)) %>%
		group_by(GenContext.norep, DomArch.norep, Lineage.final) %>%
		summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
		arrange(desc(count))
}
summ.GC.byLin <- function(x) { x %>%
		filter(!grepl("^-$", GenContext.norep)) %>%
		filter(!grepl("^-$", DomArch.norep)) %>%
		filter(!grepl("^-$", Lineage.final)) %>% filter(!grepl("^NA$", DomArch.norep)) %>%
		group_by(GenContext.norep, Lineage.final) %>% # DomArch.norep,
		summarise(count=n()) %>% # , bin=as.numeric(as.logical(n()))
		arrange(desc(count))
}
summ.GC <- function(x) { x %>%
		group_by(GenContext.norep) %>%
		summarise(totalcount=sum(count),
							totalDA=n_distinct(DomArch.norep),
							totallin=n_distinct(Lineage.final)) %>% # totallin=n_distinct(Lineage.final),
		arrange(desc(totalcount), desc(totalDA), desc(totallin)) %>%
		filter(!grepl(" \\{n\\}",GenContext.norep)) %>%
		filter(!grepl("^-$", GenContext.norep))
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

# ## counts: Function to obtain element counts (DA, GC)
# cat("Counts DAs and GCs.
# \nFor e.g.:
# query.sub$DomArch.norep %>%
#   counts(n)
# query.sub$GenContext.norep %>%
# counts(n)")

# ## elements2words: Function to break up ELEMENTS to WORDS for DA and GC
# cat("Converting DA to domains and GC to DAs.\n2 switches: da2doms and gc2da
# \nFor e.g.:
# query.sub$DA.doms <- query.sub$DomArch.norep %>%
#   elements2words(\"da2doms\")
# query.sub$GC.da <- query.sub$GenContext.norep %>%
# 	elements2words(\"gc2da\")")


# ## words2wc: Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
# cat("Word counts for broken up domains from DAs and DAs from GCs.
# \nFor e.g.:
# DA.doms.wc <- query.sub$DA.doms %>%
#   words2wc()")
