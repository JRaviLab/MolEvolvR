## Function to generate MSA: fasta —> PDF
## Created: Jun 19, 2019
## Modified: Dec 18, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse); library(here)
library(msa) # BiocManager::install("msa")
library(Biostrings)#; library(seqinr)
library(rMSA) # to implement kalign
library(pdftools); library(latexpdf); library(tools); library(tinytex) #needed?

##!! FEATURES, BUGS, NOTES thus far!!
# output generated before hand:ideally
# is determined by sliders?
# could find an output that is an image?
# look at other MSA outputs/fn / other packages to embed

#############
#### MSA ####
#############

## Sample Runs
msa_pdf(fasta_path="data/alns/pspb.gismo.fa" )#, out_path="data/msapdf")

#########################################
## Generates MSA PDF from a Fasta file ##
#########################################
msa_pdf <- function(fasta_path, out_path=NULL,
                    lowerbound=NULL, upperbound=NULL){
  #'Multiple Sequence Alignment
  #'
  #'Generates a multiple sequence alignment from a fasta file
  #'
  #'msa_pdf is a function that reads a fasta file and generates a multiple sequence alignment as
  #'a pdf
  #'
  #'@param fasta_path Character. The path location of the fasta file to be read.
  #'@param out_path Character. The path location of the output pdf to write.
  #'Default is NULL. If value is NULL, the pdf is written to the same directory as the fasta file.
  #'@param lowerbound Numeric. The column that determines the starting location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  #'@param upperbound Numeric. The column that determines the ending location of the MSA.
  #'Default is NULL. If value is NULL, the entire multiple sequence alignment is printed.
  # #'@examples msa_pdf()

  ## SAMPLE ARGUMENTS to test run
  fasta_path=here("../molevol_data/project_data/phage_defense/full_analysis_20210108/g3d.both_lin.gen.da_sub.fa")
  # out_path="../molevol_data/project_data/phage_defense/full_analysis_20210108/g3d.both_lin.gen.da.pdf"
  lowerbound=NULL; upperbound=NULL

  ## PATH DEFINITIONS
  # path+filename
  fastafile_split <- strsplit(fasta_path, "/")[[1]]
  # retrieving the filename without the '.fasta' extension
  # --> to prefix .tex and .pdf
  fastafile_name <- fastafile_split[length(fastafile_split)] %>%
    str_replace(pattern=".fasta|.fa|.faa", replacement="")
  # path to the fasta file
  inpath <- paste0(fastafile_split[1:(length(fastafile_split)-1)],
                       collapse="/")
  # tex & pdf file paths
  # tex_path <- paste0(inpath, "/", fastafile_name, ".tex")
  pdf_path <- paste0(inpath, "/", fastafile_name, ".pdf")
  print(pdf_path) ## Needs initialization
  aln_path <- paste0(inpath, "/", fastafile_name, "aln")

  ## MSA
  #!! don't want accession numbers, must be mapped to a species?
  my_seqs <- readAAStringSet(fasta_path)
  my_seqs_msa <- msa(my_seqs)

  ## Printing MSA to TeX
  #print the whole MSA if either bound is NULL
  if(is.null(lowerbound) | is.null(upperbound)){
    msaPrettyPrint(x=my_seqs_msa, output="pdf",
                   #alFile=fasta_path,
                   file=pdf_path,
                   #aesthetic
                   showNames="left", showLogo="top",
                   # showConsensus="top",
                   # showNames="none", showLogo="none",
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
                   alFile=fasta_path, file=pdf_path,
                   #if user doesn't input either lowerbound or upperbound,
                   #default to either length, or 0
                   #!!how do I get the length?
                   y=c(lowerbound,upperbound),
                   #aesthetic
                   showNames="left", showLogo="top",
                   # showConsensus="top",
                   # showNames="none", showLogo="none",
                   logoColors="rasmol",
                   shadingMode="functional",
                   shadingModeArg="structure",
                   shadingColors="blues",
                   consensusColors="ColdHot",
                   askForOverwrite=FALSE, verbose=FALSE,
                   furtherCode=c("\\defconsensus{.}{lower}{upper}",
                                 "\\showruler{1}{top}"))
  }
  system(paste("chmod g+r", fastafile_name))
  file.rename(paste0(fastafile_name, ".pdf"), pdf_path)
  system(paste("chmod g+r", pdf_path))
  ############
  ## REMOVE ##
  ############
  ## To convert TeX to PDF
  ## Errors with large number of sequences
  # print(pdf_path)
  # texi2dvi(file=tex_path, clean= TRUE, pdf=TRUE)
  # tex2pdf(tex_path, clean=TRUE) ## DOESN'T WORK EITHER!

  ### BELOW not working: maybe use function to move file after texi2pdf called?
  # curr_dir <- here("data/alns")

  ##??
  # if(is.null(out_path)){
  #   #texi2pdf outputs files to the current working directory
  #   # make it output to directory of the fasta_file
  #   # texi2pdf(file=paste0(fasta_path, "/", fastafile_name,".tex"), clean= TRUE)
  #   ##?? Alternative that works?
  #   tex2pdf(paste0(fasta_path, "/", fastafile_name,".tex"))
  # }
  # else{
  #   outfile_split=strsplit(out_path, "/")[[1]]
  #   outfile_name=outfile_split[length(outfile_split)]
  #
  #   ##!! I guess, we have these above now. please check if it works!
  #   # fasta_dir <- paste0(curr_dir,"/", fasta_path)
  #   # setwd(paste0(curr_dir,"/",output_path,"/.."))
  #
  #   #texi2pdf outputs files to the current working directory
  #   # make it output to directory of the output_path
  #   # texi2pdf(file=paste0(fasta_path, "/", fastafile_name, ".tex"), clean= TRUE)
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
  prot_aa <- readAAStringSet(path=fa_file,
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

#fasta_path <- "C:/Users/samue/Desktop/pspn.31seq.aln.txt"
#my_seqs <- readAAStringSet(fasta_path) #, format="fasta", seek.first.rec=T)
#my_seqs_msa <- msa(my_seqs)
