## code to prepare `web_blastp_hit_colnames` dataset goes here

# BLASTP and related protein BLASTs
web_blastp_hit_colnames <- c(
    "Query", "AccNum",
    "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
    "QStart", "QEnd", "SStart", "SEnd",
    "EValue", "BitScore", "PcPosOrig"
)

usethis::use_data(web_blastp_hit_colnames, overwrite = TRUE)
