# Created: April 19, 2019

library(tidyverse); library(here); library(roxygen2); library(docstring)


## Function to have Title Case
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
add_leaves <- function(input_file = "data/alignments/pspa_snf7.aln",
    lin_file = "data/pspa_snf7.txt") {
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
    #' @examples add_leaves('pspa_snf7.aln', 'pspa.txt')
    #' @details The alignment file would need two columns: 1. accession +
    #' number and 2. alignment. The protein homolog accession to lineage mapping +
    #' file should have
    #'
    #' @note Please refer to the source code if you have alternate +
    #' file formats and/or column names.

    aln <- read_tsv(input_file)
    lin <- read_tsv(lin_file)
    colnames(aln) <- c("AccNum", "Alignment")

    aln_lin <- left_join(aln, lin, by = "AccNum") %>%
        select(AccNum, Alignment,
               Species, Lineage.final)

    temp <- aln_lin %>%
        separate(Lineage.final,
                 into=c("kingdom", "phylum"),
                 sep=">", remove=F,
                 extra = "merge", fill = "left") %>%
        separate(Species, into=c("genus", "spp"),
                 sep=" ", remove=F,
                 extra = "merge", fill = "left") %>%
        mutate(leaf=paste(paste0(str_sub(temp$kingdom,
                                         start=1, end=1),
                                 str_sub(temp$phylum, 1, 6)),
                          paste0(str_sub(temp$genus, start=1, end=1),
                                 str_sub(temp$spp, start=1, end=3)),
                          temp$AccNum,
                          sep="_")) %>%
        map(leaf, to_titlecase)


    return(aln_lin)
}
