## code to prepare `ipr_cln_cols` dataset goes here

library(readr)
ipr_cln_cols <- cols(
  .default = col_character(),
  TaxID = col_double(),
  SLength = col_double(),
  StartLoc = col_double(),
  StopLoc = col_double(),
  Score = col_double(),
  Status = col_logical(),
  IPRAcc = col_logical(),
  IPRDesc = col_logical(),
  Length = col_double(),
  ID = col_logical()
)

usethis::use_data(ipr_cln_cols, overwrite = TRUE)
