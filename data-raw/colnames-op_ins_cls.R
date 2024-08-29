## code to prepare `colnames.op_ins_cls` dataset goes here
colnames.op_ins_cls <- c(
    "AccNum", "GenContext.orig",
    "DomArch.PFAM", "DomArch.orig", "DomArch.TMSIG",
    "Length", "GeneName",
    "Lineage", "Species.orig",
    "Annotation", "GI"
)

usethis::use_data(colnames.op_ins_cls, overwrite = TRUE)
