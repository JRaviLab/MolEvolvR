## code to prepare `colnames.op_ins_cls.clus2table` dataset goes here

colnames.op_ins_cls.clus2table <- c(
    "AccNum", "ClustID", "ClustName.orig",
    "GenContext.orig", "DomArch.Pfam", "DomArch.orig",
    "-", "Length", "GeneName",
    "Lineage", "Species.orig", "GCA_ID",
    "Annotation", "GI"
)

usethis::use_data(colnames.op_ins_cls.clus2table, overwrite = TRUE)
