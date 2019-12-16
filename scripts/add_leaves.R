## Functions to create leaves for MSA and Phylogenetic Tree
## Created: April 19, 2019
## Modified: Dec 12, 2019 (moved to pre-msa-tree.R)
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

#################
## Pkgs needed ##
#################
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")

###################
## Adding leaves ##
###################
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

## Function to add leaves to an alignment file
add_leaves <- function(input_file = "data/rawdata_aln/pspa_snf7.gismo.aln",
                       lin = read_tsv("data/PspA.txt"),
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
  aln <- read_tsv(aln)
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
