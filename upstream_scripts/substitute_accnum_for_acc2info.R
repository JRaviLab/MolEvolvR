#!/usr/bin/env Rscript
# this script will replace the accession numbers used in acc2info with
# the unique identifiers used for each molevolvr query

source("/data/research/jravilab/molevol_scripts/R/acc2lin.R")

args <- commandArgs(trailingOnly = TRUE)
path_acc2info <- args[1]
path_query_header_map <- args[2]
# acc2info tsv will be overwritten
path_out <- args[1]

df_acc2info <- readr::read_tsv(path_acc2info)
df_query_header_map <- readr::read_tsv(path_query_header_map)
df_acc2info_substituted <- substitute_accnum_for_acc2info(df_acc2info, df_query_header_map)
print("### df_acc2finfo_substituted")
print(df_acc2info_substituted)
readr::write_tsv(df_acc2info_substituted, file = path_out, col_names = T)
