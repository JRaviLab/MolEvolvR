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



msa_pdf <- function(fasta_filepath , lowerbound=NULL, upperbound=NULL, output_path = NULL){
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
  #'@param output_path Character. The path location of the output pdf to write.
  #'Default is NULL. If value is NULL, the pdf is written to the same directory as the fasta file.
  # #'@examples msa_pdf()

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

  fastafile_split = strsplit(fasta_filepath, "/")[[1]]
  fastafile_name = fastafile_split[length(fastafile_split)]

  #texi2pdf(file = paste0(fasta_filepath,".tex"), clean= TRUE)

  ### BELOW not working: maybe use function to move file after texi2pdf called?
  curr_dir <- getwd()

  if(is.null(output_path)){
    setwd(paste0(fasta_filepath,"/.."))
    #texi2pdf outputs files to the current working directory
    # make it output to directory of the fasta_file
    texi2pdf(file = paste0(fastafile_name,".tex"), clean= TRUE)
  }
  else{
    outfile_split = strsplit(output_path, "/")[[1]]
    outfile_name = outfile_split[length(outfile_split)]

    fasta_dir <- paste0(curr_dir,"/", fasta_filepath)

    setwd(paste0(curr_dir,"/",output_path,"/.."))


    #texi2pdf outputs files to the current working directory
    # make it output to directory of the output_path
    texi2pdf(file = paste0(fasta_dir, ".tex"), clean= TRUE)

    file.rename(paste0(fastafile_name, ".pdf"), paste0(outfile_name,".pdf"))
  }
  setwd(curr_dir)
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

