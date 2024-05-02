## code to prepare `assembly_colnames` dataset goes here

assembly_colnames <- c(
  "AssemblyID",
  "bioproject", "biosample", "wgs_master", # not used
  "RefseqCategory", "TaxID", "Spp.TaxID",
  "Species", "Spp.Strain",
  "isolate", "version_status", # not used
  "assembly_level", "release_type", # not used
  "GenomeStatus",
  "seq_rel_date", "asm_name", "submitter", # not used
  "AssemblyID.GBRS",
  "paired_asm_comp", "ftp_path", # not used
  "excluded_from_refseq", "relation_to_type_material"
) # not used

usethis::use_data(assembly_colnames, overwrite = TRUE)
