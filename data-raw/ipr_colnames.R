## code to prepare `ipr_colnames` dataset goes here

ipr_colnames <- c(
  "AccNum", "SeqMD5Digest", "SLength", "Analysis",
  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
  "Status", "RunDate", "IPRAcc", "IPRDesc"
)

usethis::use_data(ipr_colnames, overwrite = TRUE)
