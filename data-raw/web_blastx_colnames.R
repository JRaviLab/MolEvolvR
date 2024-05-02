## code to prepare `web_blastx_colnames` dataset goes here

# BLASTX
web_blastx_colnames <- c(
  "Query", "AccNum",
  "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
  "QStart", "QEnd", "SStart", "SEnd",
  "EValue", "BitScore", "PcPosOrig",
  "QSFrames"
) # specific to "blastx"

usethis::use_data(web_blastx_colnames, overwrite = TRUE)
