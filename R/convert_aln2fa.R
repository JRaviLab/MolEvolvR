# Created: April 19, 2019

library(tidyverse)
library(data.table)
#library(seqRFLP)
source("R/add_leaves.R")
# library(here)
# library(roxygen2)
# library(docstring)


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


# #takes in a datatable, not a file
# convert_to_fasta <- function(aln,output_file = NULL){
#   names <- aln$Leaf_Acc
#   sequences <- aln$Alignment
#   #library(data.table)
#   aln <- data.table(names,sequences)
#   #library(seqRFLP)
#   fasta_file <- dataframe2fas(aln, output_file)
#   return(fasta_file)
# }




# aln2tsv <- function(file_path){
#   cfile <- read_delim("data/alignments/pspc.gismo.aln", delim = " ")
#   cfile <- as.data.frame(map(cfile,function(x) gsub("\\s+", "",x)))
#   colnames(cfile) <- c("AccNum", "Alignment")
# }
