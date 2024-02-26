#!/usr/bin/env Rscript

source("/data/research/jravilab/molevol_scripts/R/fa2domain.R")
source("/data/research/jravilab/molevol_scripts/upstream_scripts/00_submit_full.R")

# parse args
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
filepath_fasta <- args[2]
blast_database <- args[3]
blast_nhits <- args[4]
blast_evalue <- args[5]
phylo <- args[6]
type <- args[7]
job_code <- args[8]
submitter_email <- args[9]
advanced_options <- args[10]

fasta <- Biostrings::readAAStringSet(filepath_fasta)

# setup paths
dir_job_results <- Sys.getenv("SLURM_SUBMIT_DIR")
dir_split_by_domain <- file.path(dir_job_results, "split_by_domain")
dir.create(dir_split_by_domain)
# change working dir to allow writing the R session in the
# split_by_domain dir in the case of failure
setwd(dir_split_by_domain)
# note: no file-extension included here since interproscan handles this
filepath_ipr_out <- "interproscan"

# 1. Exec interproscan
df_iprscan <- exec_interproscan(filepath_fasta, filepath_ipr_out)
# validate data
if (nrow(df_iprscan) < 1 | is.null(df)) {quit(save = "yes", status = 1)}

# 2. Use the results to create a domain fasta
fasta_domains <- fasta2fasta_domain(filepath_fasta, df_iprscan)
# validate data
if (!is(fasta_domains, "AAStringSet") || length(fasta_domains) < 1) {
  quit(save = "yes", status = 1)
}

# set working dir back to top level of job folder
setwd(dir_job_results)
write(fasta_domains, "domains.fa")
# just in-case, construct full path to seqs
filepath_fasta_domains <- file.path(dir_job_results, "domains.fa")

# 3. Pass domain seqs into the usual app workflow
submit_full(
  dir = dir,
  sequences = filepath_fasta_domains,
  DB = blast_database,
  NHITS = blast_nhits,
  EVAL = blast_evalue,
  phylo = phylo,
  type = type,
  job_code = job_code,
  submitter_email = submitter_email,
  advanced_options = advanced_options
)
