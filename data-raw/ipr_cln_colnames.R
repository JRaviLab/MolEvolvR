## code to prepare `ipr_cln_colnames` dataset goes here

ipr_cln_colnames <- c(
    "DB.ID", "TaxID", "AccNum.noV", "AccNum",
    "SeqMD5Digest", "SLength", "Analysis", "SignDesc",
    "StartLoc", "StopLoc", "Score", "Status", "RunDate",
    "IPRAcc", "IPRDesc", "FullAccNum", "ProteinName",
    "Length", "SourceDB", "Completeness", "Lineage",
    "Species", "Name", "ShortName", "LookupTblDesc",
    "ID", "Label"
)

usethis::use_data(ipr_cln_colnames, overwrite = TRUE)
