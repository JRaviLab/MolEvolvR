## code to prepare `domarch_lookup_colnames` dataset goes here

domarch_lookup_colnames <- c("DB.ID", "ShortName", "Description", "ID")


usethis::use_data(domarch_lookup_colnames, overwrite = TRUE)
