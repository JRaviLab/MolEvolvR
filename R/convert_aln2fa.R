# Created: April 19, 2019

# library(tidyverse)
# library(here)
# library(roxygen2)
# library(docstring)

convert_aln2fa <- function(input_file, lin_file) {
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
    #' @note Please refer to the source code if you have alternate +
    #' file formats and/or column names.

    lin <- lin_file %>% 
    select(lineage) %>%
    strsplit()
}
