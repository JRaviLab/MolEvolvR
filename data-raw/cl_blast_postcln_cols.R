## code to prepare `cl_blast_postcln_cols` dataset goes here

cl_blast_postcln_cols <- c(
    "Query", "AccNum",
    "STitle", "Species", "TaxID", "Lineage", "Lineage_long", "Lineage_long_na", "Lineage_med", "Lineage_short",
    "PcPositive", "PcIdentity", "AlnLength",
    "SAccNum", "SAllSeqID",
    "Mismatch", "GapOpen",
    "QStart", "QEnd", "QLength",
    "SStart", "SEnd", "SLength",
    "EValue", "BitScore", "PcPosOrig", "QueryName"
)

usethis::use_data(cl_blast_postcln_cols, overwrite = TRUE)
