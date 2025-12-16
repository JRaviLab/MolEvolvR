#### Utility Functions for MolEvolvR ####
## Internal helper functions used across the package.
## These functions are not exported and are intended to support
## other user-facing functions in the package.


#' Filter by Genome
#'
#' @description
#' Filters the data by genome and writes the filtered data to a file.
#'
#' @param job_dir The directory containing the job files.
#' @return None. Writes the filtered data to a file.
#'
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr arrange group_by slice filter
#' @importFrom stringr str_split
#' @importFrom tibble add_column
filterByGenome <- function(job_dir) {
  df <- read_tsv(file.path(job_dir, "cln_combined_no_dupes.tsv"))
  df <- tibble::add_column(df, Assembly = "")
  for (i in 1:nrow(df)) {
    accession <- df[i, ]$AccNum
    res <- system(paste("./find_genome.sh", accession), intern = TRUE)
    res <- str_split(res, "\t")
    tryCatch(
      {
        assembly <- res[[1]][11]
        print(assembly)
        df[i, ]$Assembly <- assembly
      },
      error = function(e) {
        print("passed")
      }
    )
  }
  df_grouped <- df %>%
    arrange(desc(PcPositive)) %>%
    group_by(Species) %>%
    slice(1)
  df <- df %>% subset(Assembly %in% df_grouped$Assembly)
  write_tsv(df, file.path(job_dir, "cln_combined_by_genome.tsv"))
}

#' Filter by Domains
#'
#' @description
#' Filters the data by domains and writes the filtered data to a file.
#'
#' @param job_dir The directory containing the job files.
#' @return None. Writes the filtered data to a file.
#'
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr arrange group_by distinct ungroup
filterByDomains <- function(job_dir) {
  df <- read_tsv(file.path(job_dir, "cln_combined_no_dupes.tsv"))
  df <- df[order(desc(df$PcPositive)), ]
  df <- df %>%
    group_by(Species) %>%
    distinct(DomArch.Pfam, .keep_all = TRUE) %>%
    ungroup()
  write_tsv(df, file.path(job_dir, "blast_combined.tsv"))
}
