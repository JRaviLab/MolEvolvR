## SCRIPT to Retrieve FASTA files w/ GenBank Protein Accession Numbers
## Created: May 13, 2019 | Janani Ravi
## Last Edited: May 14

##################
library(here)
library(tidyverse)

library(ape)
library(seqinr)
##################

sessionInfo()

## Reading input data w/ accessions
all_data <- read_tsv("data/map3773c/fur-hth.op_ins_cls.txt",
                       col_names=F)
# Filtering results by Genus
all_myco <- all_data %>%
  filter(grepl(pattern="^Mycobacterium.*", x=X10))

## READING FROM GENBANK
## Sequences from GenBank w/ accession number
seq_DNAbin <- read.GenBank("AAS06323.1")
seq_DNAbin
seq_character <- read.GenBank("AAS06323.1", as.character = TRUE)
seq_character

## Retrieving Sequences w/ List of Accession Numbers
all_myco$seq <- all_myco$X2 %>%
  as.vector() %>%
  read.GenBank()

# Checking sequence attributes
all_seq <- all_myco$seq
str(all_seq)

attributes(all_seq)
names(all_seq)
attr(all_seq, "species")

# Creating Fasta IDs
all_seq_ids <- paste(attr(all_seq, "species"),
                     names(all_seq),
                     sep="_")

all_seq_ids

# Writing Protein Fasta files by AccNums
?write.dna
write.FASTA(all_seq, file="data/map3773c/all_seqs.fa",
              #format="fasta", nbcol=6, colsep=" ", colw=10)
            append=F)
# write.dna?

## Reading and rewriting protein FASTA sequences w/ proper annotations 
all_seq_seqinr_format <- read.fasta(file = "data/map3773c/all_seqs.fa",
                                    seqtype = "AA", as.string = TRUE,
                                    forceDNAtolower = FALSE)
all_seq_seqinr_format

write.fasta(sequences = all_seq_seqinr_format,
            names = all_seq_ids,
            nbchar = 10,
            file.out = "data/map3773c/all_seq_seqinr_format.fa")

## Help
# http://www.jcsantosresearch.org/Class_2014_Spring_Comparative/pdf/week_2/Jan_13_15_2015_GenBank_part_2.pdf
