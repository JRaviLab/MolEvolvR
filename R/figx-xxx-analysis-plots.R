## Modified: Apr 8, 2019
## Created: Sep 5, 2017
## Janani Ravi (@jananiravi)
## Generic form for XXX | for now, it works with search/replace "XXX" w/ protein/domain of interest
## Adapted from https://github.com/jananiravi/psp-evolution/blob/master/scripts/201709to10-scripts/figx-xxx-analysis-plots.R
## Parsing new cleaned up op_ins_cls file
## Tested with: PspMN (2017-09), toast-rack (2017-08), https://github.com/jananiravi/psp-evolution/blob/master/scripts/201709to10-scripts/

## LIBRARIES USED
library(tidyverse)
library(UpSetR)
library(gridExtra)
source("functions-op_ins_cls.R")
source("functions-plotting.R")

## WORKING DIRECTORY
here::here()

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("as_data_frame", "tibble")

## INPUT FILES
## read in the most recent op_ins_cls file below
XXX <- read_delim("XXXxxx.txt", delim="\t",
									 escape_double=FALSE, col_names=TRUE,
									 na="NA", comment="#", trim_ws=TRUE)
XXX.DAdoms <- read_delim("Top20-domains-XXX.txt", delim="\t",
												 escape_double=FALSE, col_names=FALSE,
												 na="NA", comment="#", trim_ws=TRUE)
colnames(XXX.DAdoms) <- "domains"
XXX.DA <- read_delim("Top-XXX-domarchs.txt", delim="\t",
															escape_double=FALSE, col_names=FALSE,
															na="NA", comment="#", trim_ws=TRUE)
colnames(XXX.DA) <- "domarchs"
XXX.GCDA <- read_delim("Top-XXX-neighbors.txt", delim="\t",
														escape_double=FALSE, col_names=FALSE,
														na="NA", comment="#", trim_ws=TRUE)
colnames(XXX.GCDA) <- "domarchs"

## MAPPING CLEANED UP LINEAGES to COMBINED FILE
lineages.map <- read_delim(file="20170828.organisms2lineages.map.BAE.txt",
													 delim="\t", col_names=T, comment="#", trim_ws=T)

####################
## DATA CLEANUP ####
####################
## Adding COLUMN NAMES
colnames(XXX) <- colnames.op_ins_cls()
## GI column has given trouble in the past
XXX$GI <- as.integer(XXX$GI)

## Cleaning SPECIES and LINEAGES
XXX$Species <- XXX$Species %>%
	species.cleanup()
## Appending pre-cleaned lineage column for BAE
XXX <- left_join(XXX, lineages.map,
											 by="Species")

## Cleaning DOMAIN ARCHITECTURES ##
XXX$DomArch.norep <- XXX$SIG.TM.LADB %>%
	replace.toastrack() %>%
	remove.toastrack.ignored.doms() %>%	## Check if needed
	domarch.convert2s.forupset()	%>%			## Converts TM to TM(s) -- needed?
	repeats2s()

## Cleaning GENOMIC CONTEXTS ##
XXX$GenContext.norep <- XXX$GenContext %>%
	replace.toastrack() %>%
	repeats2s()

## Removing EMPTY ROWS from DomArch.norep and GenContext.norep
XXX.sub <- XXX %>%
	remove.empty.rows()

###############################
## COUNTS OF DAs and GCs ####
###############################
## For DOMAIN ARCHITECTURES
DomArch.ge5 <- XXX.sub$DomArch.norep %>%
	counts(5)
DomArch.all <- XXX.sub$DomArch.norep %>%
	counts(1)
## FOR GENOMIC CONTEXTS (w/o constraints on XXX DomArchs)
GenContext.ge5 <- XXX.sub$GenContext.norep %>%
	counts(5)
GenContext.all <- XXX.sub$GenContext.norep %>%
	counts(1)

#####################################
## BREAKING UP ELEMENTS -> WORDS ####
#####################################
## DOMAINS in Domain Architectures
XXX.sub$DA.doms <- XXX.sub$DomArch.norep %>%
	elements2words("da2doms")
## Domain Architectures in Genomic Contexts
XXX.sub$GC.DA <- XXX.sub$GenContext.norep %>%
	elements2words("gc2da")

###################
## WORD COUNTS ####
###################
## Counts of domains in XXX domain architectures
DA.doms.wc <- XXX.sub$DA.doms %>%
	words2wc()
# DA.doms.wc.ge5 <- filter.freq(DA.doms.wc, 5)
## Counts of domain architectures in the GenContext of XXX homologs
GC.DA.wc <- XXX.sub$GC.DA %>%
	words2wc()
# GC.DA.wc.ge50 <- filter.freq(GC.DA.wc, 50) ## Main DAs in the neighborhood
# GC.DA.wc.ge10 <- filter.freq(GC.DA.wc, 10) ## Main DAs in the neighborhood

#############################
## TOP Doms, DAs and GCs ####
#############################
## Main XXX domain architectures -- Counts by DA and Lineage
XXX.DA.summ.byLin <- summ.DA.byLin(XXX.sub)
XXX.DA.summ <- summ.DA(XXX.DA.summ.byLin)
## Main Genomic Contexts -- Summarized by XXX DomArchs & Lineage
XXX.GC.summ.byDALin <- summ.GC.byDALin(XXX.sub)
XXX.GC.summ.byLin <- summ.GC.byLin(XXX.sub)
XXX.GC.summ <- summ.GC(XXX.GC.summ.byDALin)

###########################
## Stacked Lineage Plots ##
###########################
XXX.sub <- XXX.sub %>% filter(grepl("a", Lineage.final))
# lineage.plot(XXX.sub, 10, "da2doms")
# lineage.plot(XXX.sub, 10, "gc2da")
lineage.DA.plot(XXX.sub, XXX.DA.summ.byLin) #, "da2doms"
lineage.neighbors.plot(XXX.sub, "GenContext.norep")
lineage.domain_repeats.plot(XXX.sub, "SIG.TM.LADB")

#################################
## UPSET plots for DA & Doms ####
#################################
# upset.plot.da.doms(XXX.sub, 10)
# upset.plot.gc.da(XXX.sub, 10, GenContext.norep)

upset.plot(XXX.sub, 10, "da2doms")
# upset.plot(XXX.sub, 50, "gc2da")


############################
## WRITING OUTPUT FILES ####
############################
## Last written on Sep 4, 2017 ##
# write_delim(x=XXX.sub, "XXX.sub.v2.txt",
# 						delim="\t", col_names=TRUE)

## Last written on Aug 15, 2017 ##
# write_delim(x=XXX.sub, "XXX.sub.v1.txt",
# 						delim="\t", col_names=TRUE)
# write_delim(x=DA.doms.wc, "XXX.queryDA.domains.wordcounts.txt",
# 						delim="\t", col_names=TRUE)
# write_delim(x=GC.DA.wc, "XXX.GC.DA.wordcounts.txt",
# 						delim="\t", col_names=TRUE)
# write_delim(x=DomArch.all, "XXX.queryDA.wordcounts.txt",
# 						delim="\t", col_names=TRUE)
# write_delim(x=GenContext.all, "XXX.GC.wordcounts.txt",
# 						delim="\t", col_names=TRUE)
## Last written on Aug 16
# write_delim(x=XXX.DA.summ, "XXX.DA.summ.txt",
# 						delim="\t", col_names=TRUE)
# write_delim(x=XXX.GC.summ, "XXX.GC.summ.txt",
# 						delim="\t", col_names=TRUE)
