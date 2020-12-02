## Formatting blastclust output -- create cluster identifiers based on order they appear in blastclust,
##  cluster number, and number of proteins in the cluster

## Created on 2020.05.29
## Last edit on 2020.06.17

# Load tidyverse library
library('tidyverse')
## HPC reading in arguments
args <- commandArgs(trailingOnly = TRUE)

## Local
# args <- c("data/molevolvr_outputs/WP_001901328.1_Vibrio_cholerae_out/WP_001901328.1_Vibrio_cholerae.bclust.L60S80.out",
#           "data/molevolvr_outputs/WP_001901328.1_Vibrio_cholerae_out/WP_001901328.1_Vibrio_cholerae.refseq.1e-5.cln.txt")

clust2tbl <- function(clust, blast) {
  # clust <- 'data/molevolvr_outputs/phage_defense/WP_001901328.1_Vibrio_cholerae_out/WP_001901328.1_Vibrio_cholerae.bclust.L60S80.out'
  # blast <- 'data/molevolvr_outputs/phage_defense/WP_001901328.1_Vibrio_cholerae_out/WP_001901328.1_Vibrio_cholerae.refseq.1e-5.cln.txt'
  clust_out <- read_tsv(file = clust, col_names = F)
  blast_out <- read_tsv(file = blast, col_names = T)
  ## Count the number of accession numbers in a cluster\
  # Counting number of spaces between acc. no. +1
  clust_out$NumAccs <- map(.x = clust_out$X1, function(x) { (str_count(string = x, pattern = " ")+1) })
  ## Create empty vectors to store information
  empty_vec <- c("ClusterID")
  empty_vec2 <- c("RowNum")
  ## Adding empty vectors to dataframe
  clust_out[ , empty_vec] <- NA
  clust_out[ , empty_vec2] <- NA
  ## Counting number of rows to add to the RowNum column -- used for creating cluster name
  rows <- as.numeric(rownames(clust_out)) %>%
    str_pad(width = 4, pad = 0)
  ## Add row number info to dataframe
  clust_out[ , "RowNum"] <- rows
  # Name columns
  colnam <- list('AccNum', 'NoAccs', 'ClusterID', 'NumRows')
  ## Create cluster name from row number (num of clusters) and num of accessions
  clust_out <- clust_out %>%
    `colnames<-`(colnam) %>%
    mutate(ClusterID = paste0(NumRows, ".", NoAccs))

  # Store data frame column as vector
  myvar <- c("AccNum", "ClusterID")
  ## Add only the 2 colummns wanted to a new varible
  clusters <-  clust_out[myvar]
  # Initialize empty data frame
  new_clust <- data.frame(ClusterID = character(0), AccNum = character(0), stringsAsFactors = F)

  ## Assigning each sseqid to a ClusterID
  # Iterate over dataframe
  for (i in 1:nrow(clusters)) {
    # Extract the cluster name for this row
    cname = clusters$ClusterID[i]
    # Split accession numbers by a space
    # vals is a vecto of all accession numbers in a row
    vals = clusters$AccNum[i] %>% strsplit(split = " ") %>% unlist()

    # Iterate over each element in vals
    for (v in vals) {
      # add each accession num w/ corresponding cluster ID to new df
      new_clust[nrow(new_clust)+1,] = c(cname, v) }
  }

  blast_clustnames <- merge(blast_out, new_clust, by = "AccNum")

  for (i in 1:nrow(blast_clustnames)) {
    blast_clustnames$ClusterID[i] <- str_c(blast_clustnames$ClusterID[i])
  }

  first_prot <- as.data.frame(word(clust_out$AccNum), word(clust_out$ClusterID))

  ## write the new file as a TSV file
  newarg <- gsub('.bclust.L[0-9][0-9]S[0-9][0-9].out', '', clust)

  write_tsv(new_clust, file = paste0(newarg, '.clusterIDs'), append = F)
  write_tsv(first_prot, file = paste0(newarg, ".cluster_reps"), col_names = F, append = F)
  write_tsv(blast_clustnames, file = blast, col_names = T, append = F)
}

clust2tbl(args[1], args[2])
