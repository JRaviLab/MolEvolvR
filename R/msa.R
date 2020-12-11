## Function to generate MSA: fasta —> PDF
## Created: Jun 19, 2019
## Modified: Dec 18, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
library(msa); library(Biostrings)#; library(seqinr)
library(rMSA) # to implement kalign
library(pdftools); library(latexpdf); library(tools); library(tinytex) #needed?

##!! FEATURES, BUGS, NOTES thus far!!
#output generated before hand:ideally
#is determined by sliders
#for psp, don't need to get too complex, tree
#and alignment generated backend
#could find an output that is an image?
#look at other outputs and functions of msa that produce other outputs
#or other packages, to embed

#############
#### MSA ####
#############

msa_pdf(fasta_filepath = "data/alns/pspb.gismo.fa" )#, out_filepath = "data/msapdf")
msa_pdf("data/alns/pspc.gismo.fa")


####################
msa_pdf <- function(fasta_filepath, out_filepath = NULL,
                    lowerbound=NULL, upperbound=NULL){
  #'Multiple Sequence Alignment
  #'
  #'Generates a multiple sequence alignment from a fasta file
  #'
  #'msa_pdf is a function that reads a fasta file and generates a multiple sequence alignment as
  #'a pdf
  #'
  #'@param fasta_filepath Character. The path location of the fasta file to be read.
  #'@param out_filepath Character. The path location of the output pdf to write.
  #'Default is NULL. If value is NULL, the pdf is written to the same directory as the fasta file.
  #'@param lowerbound Numeric. The column that determines the starting location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  #'@param upperbound Numeric. The column that determines the ending location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  # #'@examples msa_pdf()

  ## SAMPLE ARGUMENTS to test run
  # fasta_filepath=paste0(here("data/alns/"), dropdown_var, ".fasta")
  # out_filepath=NULL; lowerbound=NULL; upperbound=NULL

  ## PATH DEFINITIONS
  # path+filename
  fastafile_split = strsplit(fasta_filepath, "/")[[1]]
  # retrieving the filename without the '.fasta' extension
  # --> to prefix .tex and .pdf
  fastafile_name = fastafile_split[length(fastafile_split)] %>%
    str_replace(pattern=".fasta", replacement="")
  # path to the fasta file
  fasta_path <- paste0(fastafile_split[1:(length(fastafile_split)-1)],
                       collapse="/")
  # tex & pdf filepaths
  # tex_filepath <- paste0(fasta_path, "/", fastafile_name, ".tex")
  pdf_filepath <- paste0(fasta_path, "/", fastafile_name, ".pdf")
  print(pdf_filepath) ## Needs initialization

  ## MSA
  #!! don't want accession numbers, must be mapped to a species?
  my_seqs <- readAAStringSet(fasta_filepath)
  my_seqs_msa <- msa(my_seqs)

  ## Printing MSA to TeX
  #print the whole MSA if either bound is NULL
  if(is.null(lowerbound) | is.null(upperbound)){
    msaPrettyPrint(x=my_seqs_msa, output="pdf",
                   alFile=fasta_filepath, file=pdf_filepath,
                   #aesthetic
                   showNames="left", showLogo="top",
                   # showConsensus="top",
                   # showNames = "none", showLogo = "none",
                   logoColors="rasmol",
                   shadingMode="functional",
                   shadingModeArg="structure",
                   shadingColors="blues",
                   consensusColors="ColdHot",
                   askForOverwrite=FALSE, verbose=FALSE,
                   furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                 "\\showruler{1}{top}"))
  } else{ #If bounds are not NULL
    msaPrettyPrint(my_seqs_msa, output="pdf",
                   alFile=fasta_filepath, file=pdf_filepath,
                   #if user doesn't input either lowerbound or upperbound,
                   #default to either length, or 0
                   #!!how do I get the length?
                   y=c(lowerbound,upperbound),
                   #aesthetic
                   showNames="left", showLogo="top",
                   # showConsensus="top",
                   # showNames = "none", showLogo = "none",
                   logoColors="rasmol",
                   shadingMode="functional",
                   shadingModeArg="structure",
                   shadingColors="blues",
                   consensusColors="ColdHot",
                   askForOverwrite=FALSE, verbose=FALSE,
                   furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                 "\\showruler{1}{top}"))
  }

  file.rename(paste0(fastafile_name, ".pdf"), pdf_filepath)

  ############
  ## REMOVE ##
  ############
  ## To convert TeX to PDF
  ## Errors with large number of sequences
  # print(pdf_filepath)
  # texi2dvi(file=tex_filepath, clean= TRUE, pdf=TRUE)
  # tex2pdf(tex_filepath, clean=TRUE) ## DOESN'T WORK EITHER!

  ### BELOW not working: maybe use function to move file after texi2pdf called?
  # curr_dir <- here("data/alns")

  ##??
  # if(is.null(out_filepath)){
  #   #texi2pdf outputs files to the current working directory
  #   # make it output to directory of the fasta_file
  #   # texi2pdf(file = paste0(fasta_path, "/", fastafile_name,".tex"), clean= TRUE)
  #   ##?? Alternative that works?
  #   tex2pdf(paste0(fasta_path, "/", fastafile_name,".tex"))
  # }
  # else{
  #   outfile_split = strsplit(out_filepath, "/")[[1]]
  #   outfile_name = outfile_split[length(outfile_split)]
  #
  #   ##!! I guess, we have these above now. please check if it works!
  #   # fasta_dir <- paste0(curr_dir,"/", fasta_filepath)
  #   # setwd(paste0(curr_dir,"/",output_path,"/.."))
  #
  #   #texi2pdf outputs files to the current working directory
  #   # make it output to directory of the output_path
  #   # texi2pdf(file = paste0(fasta_path, "/", fastafile_name, ".tex"), clean= TRUE)
  #   ##?? Alternative that works?
  #   tex2pdf(paste0(fasta_path, "/", fastafile_name,".tex"))
  #   file.rename(paste0(fastafile_name, ".pdf"), paste0(outfile_name,".pdf"))
  #   ##!! add a small clause to make sure that .pdf.pdf get replaced with .pdf!
  # }
  # ##!! is this still needed?
  # # setwd(curr_dir)
}

####################
## Work-in-progress
## Function to generate MSA using kalign
## ref: https://rdrr.io/github/mhahsler/rMSA/man/kalign.html
## https://github.com/mhahsler/rMSA
generate_msa <- function(fa_file="", outfile=""){
  prot_aa <- readAAStringSet(filepath=fa_file,
                             format="fasta")
  prot_aa

  ## Install kalign ?rMSA_INSTALL
  ## Messed up! Reimplement from kalign.R
  ## https://github.com/mhahsler/rMSA/blob/master/R/kalign.R

  source("scripts/c2r.R")

  ## align the sequences
  al <- kalign(prot_aa) #!! won't work!
  al
}

############################
## Input files: Fasta format
#my_seqs_file <- read_tsv("data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt")
# colnames(my_seqs_file) <- c("name", "sequence")
# my_alignment <- as(my_seqs_file, "BStringSet")

#for functions, take in a file path? or maybe a "StringSet"
#Parameters: filepath? Upperbound, Lowerbound, pdfortex?

#fasta_filepath <- "C:/Users/samue/Desktop/pspn.31seq.aln.txt"
#my_seqs <- readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
#my_seqs_msa <- msa(my_seqs)

