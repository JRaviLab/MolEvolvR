## code to prepare `combo_colnames` dataset goes here

combo_colnames <- c(
    "Query", "UID", "AccNum", "Species", "TaxID", "Lineage",
    "PcPositive", "ClusterID", "QueryName",
    # "AssemblyID", "GeneName", "Description", # MISSING NOW!?!
    "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
    "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP",
    "DomArch.SMART", "DomArch.TIGR"
)

usethis::use_data(combo_colnames, overwrite = TRUE)
