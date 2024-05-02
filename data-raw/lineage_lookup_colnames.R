## code to prepare `lineage_lookup_colnames` dataset goes here

lineage_lookup_colnames <- c("TaxID", "Species", "Lineage_long", "Lineage_long_na", "Lineage_med", "Lineage_short", "Lineage")

usethis::use_data(lineage_lookup_colnames, overwrite = TRUE)
