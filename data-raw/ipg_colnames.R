## code to prepare `ipg_colnames` dataset goes here

ipg_colnames <- c(
    "IPG.ID", "Source", "NucAccNum",
    "NucStart", "NucStop", "Strand",
    "AccNum", "Description",
    "Species", "Spp.Strain", "AssemblyID"
)

usethis::use_data(ipg_colnames, overwrite = TRUE)
