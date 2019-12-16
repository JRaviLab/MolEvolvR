## Pre-requisites to generate MSA and Phylogenetic Tree
## Includes the following functions:
## convert_aln2fa, to_titlecase, add_leaves, convert_aln2tsv
## Created from add_leaves.R, convert_aln2fa.R
## Modified: Dec 12, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
library(data.table)
#library(seqRFLP)
conflicted::conflict_prefer("filter", "dplyr")

##############################
## Pre-requistive functions ##
##############################
## Function to convert alignment 'aln' to fasta format for MSA + Tree
convert_aln2fa <- function(input_file = "data/alignments/pspa_snf7.aln",
                           lin = read_tsv("data/pspa_snf7.txt"),
                           #lin_file = "data/pspa_snf7.txt",
                           output_path = NULL,
                           reduced = FALSE) {
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to an alignment file w/ accessions
  #' Genomic Contexts vs Domain Architectures.
  #' @param input_file Character. Path to file. Input tab-delimited file +
  #'  alignment file accnum & alignment.
  #'  Default is 'pspa_snf7.aln'
  #' @param lin_file Character. Path to file. Protein file with accession +
  #' number to lineage mapping.
  #' Default is 'pspa.txt'
  #' @param output_path Character. Path to the writen fasta file.
  #' Default is 'NULL'
  #' @param reduced Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
  #' Default is 'FALSE'
  #' @examples add_leaves('pspa_snf7.aln', 'pspa.txt')
  #' @details The alignment file would need two columns: 1. accession +
  #' number and 2. alignment. The protein homolog accession to lineage mapping +
  #' file should have
  #' @note Please refer to the source code if you have alternate +
  #' file formats and/or column names.

  aln <- add_leaves(input_file=input_file,
                    lin = lin,
                    #lin_file=lin_file,
                    reduced=reduced)
  names <- aln$Leaf_Acc
  sequences <- aln$Alignment
  aln <- data.table(names,sequences)

  fasta <- ""
  for(i in 1:length(aln$names)){
    fasta <- paste0(fasta,">",aln$names[i],"\n",aln$sequences[i],"\n")
  }

  if(!is.null(output_path)){
    write_file(fasta, output_path)
  }

  #fasta_file <- dataframe2fas(aln, output_path)
  return(fasta)
}

##############################
## NEEDS FIXING!
# convert_aln2tsv <- function(file_path){
#   cfile <- read_delim("data/alignments/pspc.gismo.aln", delim = " ")
#   cfile <- as.data.frame(map(cfile,function(x) gsub("\\s+", "",x)))
#   colnames(cfile) <- c("AccNum", "Alignment")
# }

##############################
## Function to convert to 'Title Case'
to_titlecase <- function(x, y=" ") {
  #' Changing case to 'Title Case'
  #' @author Andrie, Janani Ravi
  #' @description Translate string to Title Case w/ delimitter.
  #' @aliases totitle, to_title
  #' @usage to_titlecase(text, delimitter)
  #' @param x Character vector.
  #' @param y Delimitter. Default is space (" ").
  #' @seealso chartr, toupper, and tolower.
  s <- strsplit(x, y)[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=y)
}

##############################
## Function to add leaves to an alignment file
##!! Add DA to leaves?
add_leaves <- function(input_file = "data/rawdata_aln/pspa_snf7.gismo.aln",
                       lin = read_tsv("data/rawdata_tsv/all_raw.txt"), #!! finally change to all_clean.txt!!
                       #lin_file = "data/rawdata_tsv/PspA.txt",
                       reduced = FALSE) {
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to an alignment file w/ accessions
  #' Genomic Contexts vs Domain Architectures.
  #' @param input_file Character. Path to file. Input tab-delimited file +
  #'  alignment file accnum & alignment.
  #'  Default is 'pspa_snf7.aln'
  #' @param lin_file Character. Path to file. Protein file with accession +
  #' number to lineage mapping.
  #' Default is 'pspa.txt'
  #' @param reduced Boolean. If TRUE, a reduced data frame will be generated with
  #' only one sequence per lineage
  #' @examples add_leaves('pspa_snf7.aln', 'pspa.txt')
  #' @details The alignment file would need two columns: 1. accession +
  #' number and 2. alignment. The protein homolog accession to lineage mapping +
  #' file should have
  #'
  #' @note Please refer to the source code if you have alternate +
  #' file formats and/or column names.


  #1. Read aln files w/ read_file
  #2. paste and collapse files so they can be read w/ tsv
  #3. If the file has 1 column, separate it
  aln <- read_file(input_file)
  aln <- paste(aln,sep = "\\s+",collapse = "\\t")
  aln <- read_tsv(aln, col_names=F)
  if(length(aln)==1){
    colnames(aln) <- "x1"
    aln <- separate(aln, x1,c("x1","x2"),sep="\\s+")
  }

  #lin <- read_tsv(lin_file)

  colnames(aln) <- c("AccNum", "Alignment")

  aln_lin <- left_join(aln, lin, by = "AccNum") %>%
    select(AccNum, Alignment,
           Species, Lineage)

  #Removes rows with NA
  aln_lin <- aln_lin[complete.cases(aln_lin),]

  if(reduced){
    #Removes duplicate lineages
    aln_lin <- aln_lin %>% distinct(Lineage, .keep_all = TRUE)
  }

  temp <- aln_lin %>%
    separate(Lineage,
             into=c("kingdom", "phylum"),
             sep=">", remove=F,
             extra = "merge", fill = "left") %>%
    separate(Species, into=c("genus", "spp"),
             sep=" ", remove=F,
             extra = "merge", fill = "left") %>%
    #Can remove the temp and should work
    mutate(leaf=paste(paste0(str_sub(kingdom,
                                     start=1, end=1),
                             str_sub(phylum, 1, 6)),
                      paste0(str_sub(genus, start=1, end=1),
                             str_sub(spp, start=1, end=3)),
                      AccNum,
                      sep="_"))
  temp$leaf <- map(temp$leaf, to_titlecase)
  temp <- temp %>%
    mutate(Leaf_Acc = (paste(leaf,AccNum,sep="_")))

  #Combine and run through add leaves
  #3 columns ACc Sequence Leaf resutl
  #Create Leaf_AccNum pasted together
  #2 columns Leaf)AccNum  and Sequence Far left
  leaf_aln <- temp %>%
    select(Leaf_Acc,Alignment)
  return(leaf_aln)
}

##############################
## convert_accnum2fa


##############################