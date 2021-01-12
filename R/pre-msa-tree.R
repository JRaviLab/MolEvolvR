## Pre-requisites to generate MSA and Phylogenetic Tree
## Includes the following functions:
## convert_aln2fa, to_titlecase, add_leaves
## generate_all_aln2fa
## convert_aln2tsv??, convert_accnum2fa??
## Created from add_leaves.R, convert_aln2fa.R, all_aln2fa.R
## Modified: Dec 24, 2019 | Jan 2021
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(here); library(tidyverse); library(data.table)
library(rentrez)
library(msa)
library(furrr); library(future); library(doParallel)
registerDoParallel(cores=detectCores()-1)
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
                                     start=1, end=1),
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


add_name = function(data,
                    accnum_col="AccNum", spec_col="Species", lin_col="Lineage",
                    lin_sep=">", out_col="Name")
{
  #' @author Samuel Chen, Janani Ravi
  #' @description This function adds a new 'Name' column that is comprised of components from
  #' Kingdom, Phylum, Genus, and species, as well as the accession
  #' @param data Data to add name column to
  #' @param accnum_col Column containing accession numbers
  #' @param spec_col Column containing species
  #' @param lin_col Column containing lineage
  #' @param lin_sep Character separating lineage levels
  #' @param out_col Column that contains the new 'Name' derived from Species,
  #'  Lineage, and AccNum info
  #' @return Original data with a 'Name' column

  cols = c(accnum_col, "Kingdom","Phylum","Genus", "Spp")
  split_data = data %>% separate(col = lin_col, into = c("Kingdom", "Phylum"), sep = lin_sep) %>%
    mutate(Kingdom = strtrim(Kingdom,1),
           Phylum = strtrim(Phylum,6))
  if(!is.null(spec_col))
  {
    split_data =  split_data %>%
      separate(spec_col, into = c("Genus", "Spp"),sep=" ") %>%
      mutate(Genus = strtrim(Genus,1),
             Spp = word(string = Spp, start = 1))
  }
  else
  {
    split_data$Genus = ""
    split_data$Spp = ""
  }

  split_data = split_data %>%
    select(all_of(cols))
  split_data[is.na(split_data)] = ""

  accnum_sym = sym(accnum_col)

  Leaf = split_data %>%
    mutate(Leaf = paste0(Kingdom, Phylum, "_",
                         Genus, Spp, "_",
                         {{accnum_sym}} ))%>%
    pull(Leaf) %>%
    stringi::stri_replace_all_regex(pattern = "^_|_$", replacement = "")%>%
    stringi::stri_replace_all_regex(pattern = "_+", replacement = "_")


  # out_col = sym(out_col)
  data =  mutate(data, l = Leaf) %>%
    setnames(old = "l",new = out_col)
  return(data)
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


acc2fa <- function(accNum_vec, out_path, plan = "sequential")
{
  #' acc2fa converts protein accession numbers to a fasta format.
  #' Resulting fasta file is written to the out_path.
  #' @author Samuel Chen, Janani Ravi
  #' @keywords accnum, fasta
  #' @param accNum_vec Character vector containing protein accession numbers to generate fasta sequences for.
  #' Function may not work for vectors of length > 10,000
  #' @param out_path String. Location where fasta file should be written to.
  #' @example acc2fa(accNum_vec = c("ACU53894.1", "APJ14606.1", "ABK37082.1"), out_path = "my_proteins.fasta")

  if(length(accNum_vec) > 0){

    partition <- function(v, groups){
      # Partition data to limit number of queries per second for rentrez fetch:
      # limit of 10/second w/ key
      l <- length(v)

      partitioned <- list()
      for(i in 1:groups)
      {
        partitioned[[i]] <- v[seq.int(i,l,groups)]
      }

      return(partitioned)
    }

    plan(strategy = plan, .skip = T)

    min_groups = length(accNum_vec)/200
    groups <- min(max(min_groups,15) ,length(accNum_vec))

    partitioned_acc <- partition(accNum_vec, groups )
    sink(out_path)


    # a <- foreach::foreach(x = 1:length(partitioned_acc), .inorder=TRUE, .packages = "rentrez") %dopar% {
    #   cat(
    #     entrez_fetch(id = partitioned_acc[[x]],
    #                  db = "protein",
    #                  rettype = "fasta",
    #                  api_key = "YOUR_KEY_HERE"
    #     )
    #   )
    # }
    a = map(1:length(partitioned_acc), function(x)
    {
      if(x %% 9 == 0)
      {
        Sys.sleep(1)
      }
      f = future(
        entrez_fetch(id = partitioned_acc[[x]],
                     db = "protein",
                     rettype = "fasta",
                     api_key = "YOUR_KEY_HERE"
        )
      ) #%...>% cat()
    }
    )
    for(f in a)
    {
      cat(value(f))
    }

    # a <- future_map(1:length(partitioned_acc), function(x)
    # {
    #   if(x%%5 == 0)
    #   {
    #     Sys.sleep(1)
    #   }
    #   cat(
    #     entrez_fetch(id = partitioned_acc[[x]],
    #                  db = "protein",
    #                  rettype = "fasta",
    #                  api_key = "YOUR_KEY_HERE"
    #     )
    #   )
    # }
    # )
    sink(NULL)
  }
}


RepresentativeAccNums <- function(prot_data,
                                  reduced = "Lineage" ,
                                  accnum_col = "AccNum")
{
  #' Function to generate a vector of one Accession number per distinct observation from 'reduced' column
  #' @author Samuel Chen, Janani Ravi
  #' @param prot_data Data frame containing Accession Numbers
  #' @param reduced Column from prot_data from which distinct observations
  #' will be generated from.
  #' One accession number will be assigned for each of these observations
  #' @param accnum_col Column from prot_data that contains Accession Numbers

  # Get Unique reduced column and then bind the AccNums back to get one AccNum per reduced column
  reduced_sym <- sym(reduced); accnum_sym = sym(accnum_col)

  distinct_reduced <- prot_data %>% select({{reduced_sym}}) %>% distinct()

  AccNum_vect = c()

  for(i in 1:nrow(distinct_reduced))
  {
    r = toString(distinct_reduced[i,reduced])

    AccNum = toString((prot_data %>%
                         filter({{reduced_sym}} == r))[1,accnum_col])

    AccNum_vect <- append(AccNum_vect, AccNum)


  }
  return(AccNum_vect)
}

alignFasta <- function(fasta_file, tool = "Muscle", outpath = NULL)
{
  #' Perform a Multiple Sequence Alignment on a FASTA file.
  #' @author Samuel Chen, Janani Ravi
  #' @param fasta_file Path to the FASTA file to be aligned
  #' @param tool Type of alignment tool to use. One of three options: "Muscle", "ClustalOmega", or "ClustalW"
  #' @param outpath Path to write the resulting alignment to as a FASTA file. If NULL, no file is written
  #' @return aligned fasta sequence as a MsaAAMultipleAlignment object

  fasta <- readAAStringSet(fasta_file)

  aligned <- switch(
    tool,
    "Muscle" = msaMuscle(fasta),
    "ClustalOmega" = msaClustalOmega(fasta),
    "ClustalW" = msaClustalW(fasta)
  )

  if(typeof(outpath) == "character")
  {
    write.MsaAAMultipleAlignment(aligned, outpath)
  }
  return(aligned)
}

write.MsaAAMultipleAlignment <- function(alignment, outpath)
{
  #' Write MsaAAMultpleAlignment Objects as algined fasta sequence
  #' MsaAAMultipleAlignment Objects are generated from calls to msaClustalOmega
  #' and msaMuscle from the 'msa' package
  #' @author Samuel Chen, Janani Ravi
  #' @param alignment MsaAAMultipleAlignment object to be written as a fasta
  #' @param outpath Where the resulting FASTA file should be written to

  l <- length(rownames(alignment))
  fasta <- ""
  for(i in 1:l)
  {
    fasta = paste0(fasta, paste(">", rownames(alignment)[i]), "\n")
    seq <- toString(unmasked(alignment)[[i]])
    fasta <- paste0(fasta, seq, "\n")
  }
  write_file(fasta, outpath)
  return(fasta)
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
