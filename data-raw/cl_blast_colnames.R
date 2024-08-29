## code to prepare `cl_blast_colnames` dataset goes here

# pre-cleanup
cl_blast_colnames <- c(
    "Query", "SAccNum", "AccNum",
    "SAllSeqID", "STitle", "Species", "TaxID",
    "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
    "QStart", "QEnd", "QLength",
    "SStart", "SEnd", "SLength",
    "EValue", "BitScore", "PcPosOrig"
)


usethis::use_data(cl_blast_colnames, overwrite = TRUE)
