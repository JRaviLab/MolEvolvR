## code to prepare `assembly_sub_colnames` dataset goes here

assembly_sub_colnames <- c(
    "TaxID", "Spp.TaxID", "Species", "Spp.Strain",
    "RefseqCategory", "GenomeStatus",
    "AssemblyID", "AssemblyID.GBRS"
)

usethis::use_data(assembly_sub_colnames, overwrite = TRUE)
