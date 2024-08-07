## code to prepare `lineage_map_cols` dataset goes here

library(readr)
lineage_map_cols <- c(
    "double",
    "character",
    "character", "character", "character", "character", "character"
)

usethis::use_data(lineage_map_cols, overwrite = TRUE)
