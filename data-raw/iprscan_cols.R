## code to prepare `iprscan_cols` dataset goes here
library(readr)

iprscan_cols <- cols(
  .default = col_character(),
  TaxID = col_double(),
  SLength = col_double(),
  SignDesc = col_character(),
  StartLoc = col_double(),
  StopLoc = col_double(),
  Score = col_double(),
  Status = col_logical(),
  IPRAcc = col_character(),
  IPRDesc = col_character(),
  Length = col_double(),
  ShortName = col_character(),
  LookupTblDesc = col_character(),
  ID = col_character(),
  Label = col_character()
)

usethis::use_data(iprscan_cols, overwrite = TRUE)
