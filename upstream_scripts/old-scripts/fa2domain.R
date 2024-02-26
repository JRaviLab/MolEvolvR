#!/usr/bin/env R
### begin code from fa2domain.sh
# adapted from original fa2domain.sh script
# not sure if this is needed since there are 
# two successive `cd` executions without any 
# executions that were sensitive to the working directory
# aside from the second `cd` call which could potentially
# be relative
dir_out <- Sys.getenv("SLURM_SUBMIT_DIR")
setwd(dir_out)

args <- commandArgs(trailingOnly = TRUE)
filepath_input <- args[1]
type <- args[2]
phylo <- args[3]
dir <- args[4]
prefix <- basename(filepath_input) |>
  fs::path_ext_remove()
setwd(dir)

cmd_iprscan <- stringr::str_glue(
  "/data/research/jravilab/molevol_scripts/upstream_scripts/04a_iprscan.sh ",
  "{filepath_input} {prefix} {dir}"
)


file.create(
  file.path(dir, stringr::str_glue("{prefix}.domains.fa"))
)
filepath_iprscan <- file.path(dir, stringr::str_glue("{prefix}.iprscan.tsv"))
filepath_fa_domains <- file.path(dir, stringr::str_glue("{prefix}.domains.fa"))

### end code from fa2domain.sh
### begin fa2domain.R
# converting previous args to variables available in combined script
#args[1] <- filepath_input
#args[2] <- file.path(dir, stringr::str_glue("{prefix}.iprscan.tsv"))
#args[3] <- file.path(dir, stringr::str_glue("{prefix}.domains.fa"))
#args[4] <- dir
#args[5] <- phylo
#args[6] <- type

# new
get_sequences(
  filepath_input,
  dir = dir,
  separate = FALSE,
  acc_file_path = "starting_accs.txt"
)
fasta <- readAAStringSet(filepath_input)
df_iprscan <- read.csv(filepath_iprscan, sep = "\t", header = FALSE)
filepath_fa_domains <- args[3]
for (item in names(in_fa)) {
  accession <- unlist(strsplit(item, " "))[1]
  pfam_count <- 1
  gene3d_count <- 1
  for (i in 1:nrow(in_ipr)) {
    if (in_ipr[i, 1] == accession && (in_ipr[i, 4] == "Pfam")) {
      sequence <- toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0(">", accession, "_", in_ipr[i, 4], "-", pfam_count, " ", in_ipr[i, 5], " ", in_ipr[i, 6])
      complete_seq <- paste(header, sequence, sep = "\n")
      cat(complete_seq, file = out_file, append = TRUE, sep = "\n")
      pfam_count <- pfam_count + 1
    } else if (in_ipr[i, 1] == accession && in_ipr[i, 4] == "Gene3D") {
      sequence <- toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0(">", accession, "_", in_ipr[i, 4], "-", gene3d_count, " ", in_ipr[i, 5], " ", in_ipr[i, 6])
      complete_seq <- paste(header, sequence, sep = "\n")
      cat(complete_seq, file = out_file, append = TRUE, sep = "\n")
      gene3d_count <- gene3d_count + 1
    } else {
    }
  }
}

submit_full(dir = args[4], sequences = out_file, phylo = args[5], by_domain = "TRUE", domain_starting = args[1], type = args[6])
