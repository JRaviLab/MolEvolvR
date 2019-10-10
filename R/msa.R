library(tidyverse)
library(msa)
# library(seqinr)
library(Biostrings)
library(tools)

#output generated before hand:ideally y
#is determined by sliders
#for psp, don't need to get too complex, tree
#and alignment generated backend
#could find an output that is an image?
#look at other outputs and functions of msa that produce other outputs
#or other packages, to embed



msa_pdf <- function(fasta_filepath , lowerbound=NULL, upperbound=NULL){
  #'Multiple Sequence Alignment
  #'
  #'Generates a multiple sequence alignment from a fasta file
  #'
  #'msa_pdf is a function that reads a fasta file and generates a multiple sequence alignment as
  #'a pdf
  #'
  #'@param fasta_filepath Character. The path location of the fasta file to be read.
  #'@param lowerbound Numeric. The column that determines the starting location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  #'@param upperbound Numeric. The column that determines the ending location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  # #'@examples msa_pdf()

  #only works for a fasta file
  #first convert tsv to fastafile: pretty easy but something needs to be done first
  #tab replaced with something and < in front
  # don't want accesation numbers, must be mapped to a species
  #This is
  my_seqs <- readAAStringSet(fasta_filepath)
  my_seqs_msa <- msa(my_seqs)

  #print the whole MSA if either bound is NULL
  if(is.null(lowerbound) | is.null(upperbound)){
    msaPrettyPrint(my_seqs_msa, output="tex",
                   file=paste0(fasta_filepath, ".tex"),

                   #aesthetic
                   showNames="left", showLogo="top",
                   # showNames = "none",
                   # showLogo = "none",
                   logoColors="rasmol", # “chemical”, “rasmol”, “hydropathy”, “structure”, “standard area”, “accessible area”
                   shadingMode="functional", # or "similar"
                   shadingModeArg="structure",
                   shadingColors="blues",
                   consensusColors="ColdHot",
                   askForOverwrite=FALSE, verbose=FALSE,
                   furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                 "\\showruler{1}{top}"))
  }
  #If bounds are not NULL
  else{
    msaPrettyPrint(my_seqs_msa, output="tex",
                   file=paste0(fasta_filepath, ".tex"),

                   #if user doesn't input either lowerbound or upperbound, default to either length, or 0
                   #how do I get the length?
                   y=c(lowerbound,upperbound),
                   #aesthetic
                   ##showNames="left", showLogo="top",
                   showNames = "none",
                   showLogo = "none",
                   ##logoColors="rasmol", # “chemical”, “rasmol”, “hydropathy”, “structure”, “standard area”, “accessible area”
                   shadingMode="functional", # or "similar"
                   shadingModeArg="structure",
                   shadingColors="blues",
                   consensusColors="ColdHot",
                   askForOverwrite=FALSE, verbose=FALSE,
                   furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                 "\\showruler{1}{top}"))
  }

  texi2pdf(file = paste0(fasta_filepath,".tex"), clean= TRUE)
}

## Input files: Fasta format
#my_seqs_file <- read_tsv("data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt")
# colnames(my_seqs_file) <- c("name", "sequence")
# my_alignment <- as(my_seqs_file, "BStringSet")

#for functions, take in a file path? or maybe a "StringSet"
#Parameters: filepath? Upperbound, Lowerbound, pdfortex?

#fasta_filepath <- "C:/Users/samue/Desktop/pspn.31seq.aln.txt"
#my_seqs <- readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
#my_seqs_msa <- msa(my_seqs)

