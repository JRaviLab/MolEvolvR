#!/usr/bin/env Rscript
t_0 <- Sys.time()
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

msg_args <- stringr::str_glue(
  "args from split_by_domain-runner:\n",
  "\tpaste(args, collapse = ',')\n"
)
print(msg_args)

fasta <- Biostrings::readAAStringSet(filepath_fasta)

# setup paths
dir_job_results <- Sys.getenv("SLURM_SUBMIT_DIR")

# wrap split by domain code in function
# since this entire expression needs to be wrapped in tryCatch
# (to handle the case of no domains returned from the selected analyses)
# just to minimize the amount of indentation formatting
# as an alternative to wrapping the whole expression in tryCatch directly
split_by_domain <- function() {
  dir_split_by_domain <- file.path(dir_job_results, "split_by_domain")
  dir.create(dir_split_by_domain)
  # change working dir to allow writing the R session in the
  # split_by_domain dir in the case of failure
  setwd(dir_split_by_domain)
  # note: no file-extension included here since interproscan handles this
  filepath_ipr_out <- "interproscan"

  # 1. exec interproscan
  df_iprscan <- exec_interproscan(filepath_fasta, filepath_ipr_out)
  print("### df_iprscan")
  print(df_iprscan)
  # validate data
  if (nrow(df_iprscan) < 1 || is.null(df)) {
    save.image("debug-session.rda")
    # raise error to launch regular submission
    stop("validation of df_iprscan failed")
  }

  # 2. use the results to create a domain fasta
  fasta_domains <- fasta2fasta_domain(fasta, df_iprscan, verbose = TRUE)
  print("### fasta domains")
  print(fasta_domains)
  # validate data
  if (!is(fasta_domains, "AAStringSet") || length(fasta_domains) < 1) {
    save.image("debug-session.rda")
    # raise error to launch regular submission
    stop("validation of 'fasta_domains' failed")
  }

  # construct full path to domain seqs and write data
  filepath_fasta_domains <- file.path(dir_split_by_domain, "domains.fa")
  Biostrings::writeXStringSet(fasta_domains, filepath_fasta_domains)

  # set working dir back to top level of job folder
  setwd(dir_job_results)

  # 3. append domain seqs to input fasta and submit to usual app workflow
  fa_original_txt <- readLines(filepath_fasta)
  fa_domains_txt <- readLines(filepath_fasta_domains)
  fa_concat_txt <- c(fa_original_txt, fa_domains_txt)
  filepath_fasta_concat <- file.path(
    paste0(filepath_fasta |> fs::path_ext_remove(), "-w-domain_seqs.fa")
  )
  writeLines(fa_concat_txt, con = filepath_fasta_concat)

  submit_full(
    dir = dir,
    sequences = filepath_fasta_concat,
    DB = blast_database,
    NHITS = blast_nhits,
    EVAL = blast_evalue,
    phylo = phylo,
    type = type,
    job_code = job_code,
    submitter_email = submitter_email,
    advanced_options = advanced_options
  )
}

# upon any error (such as those explicitly raised by
# lacking domain data), use the original protein data
result_split_by_domain <- tryCatch(
  # try getting domains from proteins and submit them to the app
  expr = {
    split_by_domain()
    TRUE
  },
  # upon error submit the original protein data
  error = function(e) {
    print(e)
    msg <- stringr::str_glue(
      "* An error was raised during the split by domain phase.\n",
      "* Trying to submit the original data (without processing domains)\n"
    )
    warning(msg)
    submit_full(
      dir = dir,
      sequences = filepath_fasta,
      DB = blast_database,
      NHITS = blast_nhits,
      EVAL = blast_evalue,
      phylo = phylo,
      type = type,
      job_code = job_code,
      submitter_email = submitter_email,
      advanced_options = advanced_options
    )
    FALSE
  }
)

### logging
# if the domain fasta was written, count number of domain seqs
n_domains_split <- ifelse(
  file.exists(file.path(dir_job_results, "split_by_domain", "domains.fa")),
  readAAStringSet(file.path(dir_job_results, "split_by_domain", "domains.fa")) |> length(),
  NA
)
df_log_split_by_domain <- tibble::tibble(
  # format the date times to match the other molevolvr log files
  "START_DT" = t_0 |> format("%d/%m/%Y-%H:%M:%S"),
  "STOP_DT" = Sys.time() |> format("%d/%m/%Y-%H:%M:%S"),
  "job_code" = job_code,
  "was_domain_split_successful" = result_split_by_domain,
  "n_domains_split" = n_domains_split
)
readr::write_tsv(
  df_log_split_by_domain,
  file = file.path(dir_job_results, "split_by_domain", "logfile.tsv")
)
