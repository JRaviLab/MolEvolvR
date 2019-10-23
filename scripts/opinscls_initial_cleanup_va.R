library(tidyverse)

#reading data file
pspa_data <- read_delim("all.op_ins_cls",
                        "\t", escape_double = FALSE, col_names = FALSE,
                        col_types = cols(X13 = col_skip(), X14 = col_skip(),
                                         X2 = col_skip(), X7 = col_skip()),
                        trim_ws = TRUE)
pspa_data=data.frame(pspa_data)
names(pspa_data)=c("Acc", "arch", "operons", "archpfam", "archprofdb", "len","gene.name","tax","org","GCA")


#read grep domain list
domains_forgrep <- read_csv("domains_forgrep.txt",
                            col_names = FALSE)
domains_forgrep=data.frame(domains_forgrep)
colnames(domains_forgrep) <- c("grepcat")
domains_forgrep

#getting only those in domains_forgrep
pspa_data=pspa_data[grep(pattern = paste((domains_forgrep$grepcat), collapse = "|"), pspa_data$arch),]


#sorting by architecture length
pspa_data=pspa_data[order(nchar(pspa_data$arch),decreasing = T),]

#making accessions non-redundant
pspa_data=pspa_data[!duplicated(pspa_data$Acc), ]

#splitting and adding $taxend
te=sapply(strsplit(as.character(pspa_data$tax),">"),tail,1)
pspa_data$taxend = te