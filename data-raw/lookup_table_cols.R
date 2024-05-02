## code to prepare `lookup_table_cols` dataset goes here

library(readr)

lookup_table_cols <- cols(
  DB.ID = col_character(),
  ShortName = col_character(),
  Description = col_character(),
  ID = col_character()
)

usethis::use_data(lookup_table_cols, overwrite = TRUE)
