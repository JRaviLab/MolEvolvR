## Pre-requisites to generate MSA and Phylogenetic Tree
## Includes the following functions:
## convert_aln2fa, to_titlecase, add_leaves
## generate_all_aln2fa
## convert_aln2tsv??, convert_accnum2fa??
## Created from add_leaves.R, convert_aln2fa.R, all_aln2fa.R
## Modified: Dec 24, 2019
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(here); library(tidyverse)
library(data.table)
#library(seqRFLP)
conflicted::conflict_prefer("filter", "dplyr")

##############################
## Pre-requisite functions ##
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

################################
## Function to add leaves to an alignment file
##!! Add DA to leaves?
add_leaves <- function(aln_file = "",
                       lin_file = "data/rawdata_tsv/all_semiclean.txt", #!! finally change to all_clean.txt!!
                       #lin_file = "data/rawdata_tsv/PspA.txt",
                       reduced = FALSE) {
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to an alignment file w/ accessions
  #' Genomic Contexts vs Domain Architectures.
  #' @param aln_file Character. Path to file. Input tab-delimited file +
  #'  alignment file accnum & alignment.
  #'  Default is 'pspa_snf7.aln'
  #' @param lin_file Character. Path to file. Protein file with accession +
  #' number to lineage mapping.
  #' Default is 'pspa.txt'
  #' @param reduced Boolean. If TRUE, a reduced data frame will be generated with
  #' only one sequence per lineage. Default is FALSE.
  #' @examples add_leaves('pspa_snf7.aln', 'pspa.txt')
  #' @details The alignment file would need two columns: 1. accession +
  #' number and 2. alignment. The protein homolog accession to lineage mapping +
  #' file should have
  #'
  #' @note Please refer to the source code if you have alternate +
  #' file formats and/or column names.

  ## SAMPLE ARGS
  # aln_file <- "data/rawdata_aln/pspc.gismo.aln"
  # lin_file <- "data/rawdata_tsv/all_semiclean.txt"
  # reduced=F

  #1. Read aln & lineage master files files w/ read_file/read_tsv
  #2. paste and collapse files so they can be read w/ tsv
  #3. If the file has 1 column, separate it
  aln <- read_file(aln_file)
  lin <- read_tsv(lin_file)
  aln <- paste(aln,sep = "\\s+",collapse = "\\t")
  aln <- read_tsv(aln, col_names=F)
  if(length(aln)==1){
    colnames(aln) <- "x1"
    aln <- separate(aln, col=x1,
                    into=c("x1","x2"),
                    sep="\\s+")
  }

  colnames(aln) <- c("AccNum", "Alignment")

  aln_lin <- left_join(aln, lin, by = "AccNum") %>%
    select(AccNum, Alignment,
           Species, Lineage)

  #Removes rows with NA
  aln_lin <- aln_lin[complete.cases(aln_lin),]
  #Removes duplicated rows
  aln_lin <- aln_lin %>% distinct()

  ##!! REVISE or REMOVE??? MSA doesn't work with too many sequences!!
  ##!! FIX ASAP.
  if(reduced){
    #Removes duplicate lineages
    aln_lin <- aln_lin %>% distinct(Lineage, .keep_all = TRUE)
  }

  temp <- aln_lin %>%
    separate(Lineage,
             into=c("Kingdom", "Phylum"),
             sep=">", remove=F,      ### Some  lineages don't have a phylum, What should those be turned into
             extra = "merge", fill = "right") %>%
    replace_na(replace = list(Kingdom = "", Phylum = "")) %>%
    separate(Species, into=c("Genus", "Spp"),
             sep=" ", remove=F,
             extra = "merge", fill = "left") %>%

    # 3char from kingdom, 6 char from phylum, 1 char from Genus, 3 char from species
    # kingdomPhylum_GenusSpecies
    mutate(Leaf=paste(paste0(str_sub(Kingdom,
                                     start=1, end=3),
                             str_sub(Phylum, 1, 6)),
                      paste0(str_sub(Genus, start=1, end=1),
                             str_sub(Spp, start=1, end=3)),
                      # AccNum,
                      sep="_"))
  temp$Leaf <- map(temp$Leaf, to_titlecase)
  temp <- temp %>%
    mutate(Leaf_Acc = (paste(Leaf, AccNum, sep="_")))

  #Combine and run through add leaves
  #3 columns AccNum Sequence Leaf result
  #Create Leaf_AccNum pasted together
  #2 columns Leaf_AccNum and Sequence Far left
  leaf_aln <- temp %>%
    select(Leaf_Acc, Alignment)
  return(leaf_aln)
}

################################
## Function to convert alignment 'aln' to fasta format for MSA + Tree
convert_aln2fa <- function(aln_file = "",
                           lin_file = "data/rawdata_tsv/all_semiclean.txt", #!! finally change to all_clean.txt!!
                           fa_outpath = "",
                           reduced = FALSE) {
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to an alignment file w/ accessions
  #' Genomic Contexts vs Domain Architectures.
  #' @param aln_file Character. Path to file. Input tab-delimited file +
  #'  alignment file accnum & alignment.
  #'  Default is 'pspa_snf7.aln'
  #' @param lin_file Character. Path to file. Protein file with accession +
  #' number to lineage mapping.
  #' Default is 'pspa.txt'
  #' @param fa_outpath Character. Path to the written fasta file.
  #' Default is 'NULL'
  #' @param reduced Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
  #' Default is 'FALSE'
  #' @examples add_leaves('pspa_snf7.aln', 'pspa.txt')
  #' @details The alignment file would need two columns: 1. accession +
  #' number and 2. alignment. The protein homolog accession to lineage mapping +
  #' file should have
  #' @note Please refer to the source code if you have alternate +
  #' file formats and/or column names.

  ## SAMPLE ARGS
  # aln_file <- "data/rawdata_aln/pspc.gismo.aln"
  # lin_file <- "data/rawdata_tsv/all_semiclean.txt"
  # reduced=F
  # fa_outpath="data/alns/pspc.fasta"

  ## Add leaves
  aln <- add_leaves(aln=aln_file,
                    lin=lin_file,
                    reduced=reduced)
  names <- aln$Leaf_Acc
  sequences <- aln$Alignment
  aln <- data.table(names, sequences)

  ## Convert to Fasta
  fasta <- ""
  for(i in 1:length(aln$names)){
    fasta <- paste0(fasta,">",aln$names[i],"\n",aln$sequences[i],"\n")
  }

  if(!is.null(fa_outpath)){
    write_file(fasta, fa_outpath)
  }

  #fasta_file <- dataframe2fas(aln, output_path)
  return(fasta)
}

################################
## generate_all_aln2fa
generate_all_aln2fa <- function(aln_path=here("data/rawdata_aln/"),
                                fa_outpath=here("data/alns/"),
                                lin_file=here("data/rawdata_tsv/all_semiclean.txt"),
                                reduced=F)
{
  #' Adding Leaves to an alignment file w/ accessions
  #' @author Janani Ravi
  #' @keywords alignment, accnum, leaves, lineage, species
  #' @description Adding Leaves to all alignment files w/ accessions & DAs?
  #' @param aln_path Character. Path to alignment files.
  #' Default is 'here("data/rawdata_aln/")'
  #' @param lin_file Character. Path to file. Master protein file with AccNum & lineages.
  #' Default is 'here("data/rawdata_tsv/all_semiclean.txt")'
  #' @param fa_outpath Character. Path to the written fasta file.
  #' Default is 'here("data/alns/")'.
  #' @param reduced Boolean. If TRUE, the fasta file will contain only one sequence per lineage.
  #' Default is 'FALSE'.
  #' @examples generate_all_aln2fa()
  #' @details The alignment files would need two columns separated by spaces: 1. AccNum and 2. alignment. The protein homolog file should have AccNum, Species, Lineages.
  #' @note Please refer to the source code if you have alternate + file formats and/or column names.

  library(here)
  library(tidyverse)
  # aln_path <- here("data/rawdata_aln/")
  # outpath <- here("data/alns/")
  # lin_file <- here("data/rawdata_tsv/all_semiclean.txt")

  aln_filenames <- list.files(path=aln_path, pattern="*.aln")
  aln_filepaths <- paste0(aln_path,"/", aln_filenames)
  variable <- str_replace_all(basename(aln_filenames),
                              pattern=".aln", replacement="")

  ## Using purrr::pmap
  aln2fa_args <- list(aln_file=aln_filepaths,
                      fa_outpath=paste0(fa_outpath, "/", variable, ".fa"))
  pmap(.l=aln2fa_args, .f=convert_aln2fa,
       lin_file=lin_file,
       reduced=reduced)
}

################################
## convert_accnum2fa
#######
## 1 ##
#######
## Accnum2fa
# ref <- c("U15717", "U15718", "U15719", "U15720",
#          "U15721", "U15722", "U15723", "U15724")
# ref_gb <- read.GenBank(ref)
# cbind(attr(ref_gb, "species"), names(ref_gb))
# attr(ref_gb, "description")

#######
## 2 ##
#######
## GenBank 2 fasta file format
# seqinr::gb2fasta(source.file="", destination.file="")

#######
## 3 ##
#######
## Ref: https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
# retrieveseqs <- function(seqnames,acnucdb)
# {
#   myseqs <- list()   # Make a list to store the sequences
#   require("seqinr")  # This function requires the SeqinR R package
#   choosebank(acnucdb)
#   for (i in 1:length(seqnames))
#   {
#     seqname <- seqnames[i]
#     print(paste("Retrieving sequence",seqname,"..."))
#     queryname <- "query2"
#     query <- paste("AC=",seqname,sep="")
#     query(`queryname`,`query`)
#     seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
#     myseqs[[i]] <- seq
#   }
#   closebank()
#   return(myseqs)
# }
# seqnames <- c("Q10572","E3M2K8","Q8WS01","E1FUV2","A8NSK3","Q9VT99")
# seqs <- retrieveseqs(seqnames,"swissprot")

################################
## convert_aln2tsv
## NEEDS FIXING!
# convert_aln2tsv <- function(file_path){
#   cfile <- read_delim("data/alignments/pspc.gismo.aln", delim = " ")
#   cfile <- as.data.frame(map(cfile,function(x) gsub("\\s+", "",x)))
#   colnames(cfile) <- c("AccNum", "Alignment")
# }
