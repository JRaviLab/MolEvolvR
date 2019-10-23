# Adding column names to our final master file!
# Created: JR, VA
# Oct 22, 2019

library(tidyverse)
pspa_data <- read_delim("data/rawdata_opinscls/all.op_ins_cls",
                        delim="\t", escape_double=FALSE,
                        col_names=FALSE, trim_ws=TRUE,
                        col_types=cols(X13=col_skip(), X14=col_skip(),
                                         X2=col_skip(), X7=col_skip()))
## skipping 4 columns
# 2nd (ClustNum),7th (Blank),13th (Annotation) and 14th (GI)
# MISSING: TaxID & Query

## Need to convert to DF structure?
# pspa_data <- data.frame(pspa_data)

## Adding col names
# names(pspa_data)=c("Acc", "arch", "operons", "archpfam", "archprofdb",
#                    "len","gene.name","tax","org","GCA")
names(pspa_data) <- c("AccNum", "ClustName", "GenContext.orig",
                   "DomArch.Pfam", "DomArch.orig",
                   "Length", "GeneName", "Lineage", "Species.orig", "GCA_ID")

## have asked V the following Qs: arch=ClustName and archprofdb=DomArch.orig?
## query deduced from ClustName based on the list below?
## TaxID for these AccNums

write_tsv(pspa_data, "data/rawdata_tsv/all_raw_notax.txt", col_names=T)

## Core Domains For grep
# DUF1700-alpha-helical
# DUF1707-SHOCT-bihelical
# LiaI-LiaF-TM
# Toast-rack
# PspA
# PspB
# PspC
# PspM
# PspN_N
# DUF3046
# Snf7