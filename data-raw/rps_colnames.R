## code to prepare `rps_colnames` dataset goes here

rps_colnames <- c(
    "AccNum", "DB.ID", "DBSeqID",
    "PcIdentity.Dom", "PcPosOrig.Dom", # "PcPos.Dom", # Ppos missing
    "AlnLength", "Mismatch",
    "SStart", "SEnd", "DStart", "DEnd",
    "EValue", "BitScore", "TaxID"
) # TaxID missing (NA); remove?

usethis::use_data(rps_colnames, overwrite = TRUE)
