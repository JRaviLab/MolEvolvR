## code to prepare `web_blastp_desc_colnames` dataset goes here

## Downloaded as Descriptions csv
# BLASTP and related protein BLASTs
web_blastp_desc_colnames <- c(
    "Description", "Species", "CommonName", "TaxID",
    "BitScore", "TotalScore",
    "PcQCover", "EValue", "PcIdentity",
    "SLen", "AccNum"
)
# Ref: https://ncbiinsights.ncbi.nlm.nih.gov/2020/11/23/blast-new-columns/
# Description,	Scientific Name,	Common Name,	Taxid,
# Max Score,	Total Score,
# Query Cover,	E value,	Per. ident,
# Acc. Len	Accession

usethis::use_data(web_blastp_desc_colnames, overwrite = TRUE)
