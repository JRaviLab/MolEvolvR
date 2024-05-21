suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
options(readr.show_col_types = FALSE, readr.show_types = FALSE) # silence read tsv col types
library(yaml)

# for now, we're using an env var, COMMON_SRC_ROOT, to specify this folder since
# the working directory is changed in many parts of the current molevolvr
# pipeline.
common_root <- Sys.getenv('COMMON_SRC_ROOT')
source(file.path(common_root, "molevol_scripts", "R", "assign_job_queue.R"))

get_sequences <- function(sequences, acc_file_path = "accs.txt", dir_path = "~", separate = TRUE) {
  seqs <- readAAStringSet(sequences)
  cln_names <- c()
  for (accnum in names(seqs)) {
    if (grepl("\\|", accnum)) {
      accnum_cln <- strsplit(accnum, "\\|")[[1]][2]
      accnum_cln <- strsplit(accnum_cln, " ")[[1]]
    } else {
      accnum_cln <- strsplit(accnum, " ")[[1]][1]
    }
    cln_names <- append(cln_names, accnum_cln)
    write(accnum_cln, file = acc_file_path, append = TRUE)
    if (separate) {
      write(paste0(dir_path, "/", accnum_cln, ".faa"), file = "input.txt", append = TRUE)
      write(paste0(">", accnum_cln), file = paste0(accnum_cln, ".faa"), append = TRUE)
      write(toString(seqs[accnum]), file = paste0(accnum_cln, ".faa"), append = TRUE)
    }
  }
  names(seqs) <- cln_names
  writeXStringSet(seqs, sequences, format = "fasta")
  return(length(seqs))
}

make_job_name <- function(job_code, suffix = "molevol_analysis") {
  if (!is.null(job_code)) paste(job_code, suffix, sep="_") else suffix
}

#' Produces arguments to make Slurm send job status emails to the submitter.
#' The result of this function is intended to be used as an argument to sbatch.
#' If the submitter's email is NULL or empty, or if get_slurm_mails is FALSE,
#' this function returns an empty string.
#' 
#' @param submitter_email The email address of the job submitter
#' @param get_slurm_mails Whether to use Slurm's built-in email notifications
#' @param mailtypes A comma-delimited list of Slurm event types for which to send notifications
#' 
#' @return A string containing the arguments to pass to sbatch to enable email notifications, if applicable
#' 
make_email_args <- function(submitter_email, get_slurm_mails, mailtypes="END,FAIL") {
  if (get_slurm_mails && !is.null(submitter_email) && stringr::str_trim(submitter_email) != "") {
    stringr::str_glue(
      "--mail-type={mailtypes} --mail-user=\"{submitter_email}\""
    )
  } else { "" }
}

#' Executes 'cmd', returning a list of the form {result, exit_code},
#' where 'result' is the output of the command and 'exit_code' is the
#' exit code of the command. If the command fails to execute, 'exit_code'
#' will be the error code of the system call.
#' 
#' See the R system() docs for more info on the result object,
#' https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
#' 
#' @param cmd A string containing the command to execute
#' @return A list containing the result of the command and the exit code
exec_cmd <- function(cmd) {
  cmd_errored <- FALSE

  result <- tryCatch(
    system(cmd, intern=TRUE),
    error = function(e) { cmd_errored <<- TRUE }
  )
  result_status <- attr(result, 'status')

  exit_code <- ifelse(
    cmd_errored,
    result_status,
    ifelse(is.null(result_status), 0, result_status)
  )

  return(list(result=result, exit_code=exit_code))
}

# print info on slurm submission errors
# parses and returns the slurm job id if successful
submit_and_log <- function(cmd, exit = FALSE) {
  result <- exec_cmd(cmd)

  cat(file=stderr(), paste0("Ran command ", cmd, ", result: ", result$exit_code, "\n"))
  flush.console()

  if (result$exit_code != 0L) {
    if (exit == TRUE) stop(paste0("failed to submit job; error code: ", result$exit_code, "\n", "cmd=", cmd))
  }
  else {
    # parse out the job id from the result.
    # the response is always of the form 'Submitted batch job 123456',
    # so we split on spaces and take the last element as an integer.
    job_id <- tail(strsplit(result$result, " ")[[1]], n=1)
    return(as.integer(job_id))
  }
}

# if submitter_email is not NULL or empty, schedule a 'summary' job to notify
# the user via email when all the job_ids have completed or failed
submit_summary_job <- function(job_ids, submitter_email, job_dir, job_code, dep_jobs) {
  if (!is.null(submitter_email) && stringr::str_trim(submitter_email) != "") {
    type <- 'notify'
    dep_jobs <- paste(job_ids, collapse=",")
    submit_and_log(stringr::str_glue(
      "sbatch --partition NotifyQ --dependency=afterany:", paste(job_ids, collapse=":"),
      " --job-name {make_job_name(job_code, type)}",
      " /data/research/jravilab/molevol_scripts/upstream_scripts/99_send_completion_mail.R",
      " '{job_dir}' '{submitter_email}' '{job_code}' '{dep_jobs}'"
    ))
  }
}

submit_full <- function(dir = "/data/scratch", DB = Sys.getenv("BLAST_DB", unset = "refseq"), NHITS = Sys.getenv("BLAST_HITS", unset = 100), EVAL = Sys.getenv("BLAST_EVALUE", unset = 0.00001), sequences = "~/test.fa", phylo = "FALSE", by_domain = "FALSE", domain_starting = "~/domain_seqs.fa", type = "full", job_code=NULL, submitter_email=NULL, advanced_options=NULL, get_slurm_mails=FALSE) {
  # submits jobs for fasta, MSA, or accession number type submissions
  setwd(dir)

  advanced_options_names <- names(advanced_options[advanced_options==TRUE])

  # write job submission params to file
  job_args <- list(
    submission_type = type,
    database = ifelse(phylo == FALSE, DB, NA),
    nhits = ifelse(phylo == FALSE, NHITS, NA),
    evalue = ifelse(phylo == FALSE, EVAL, NA),
    submitter_email = submitter_email,
    advanced_options = advanced_options_names,
    job_code = job_code
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  # keep track of the job ids so we can join on them for the summary job
  job_ids <- c()

  num_runs <- 0
  write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")

  # assign destPartition to a known partition name just in case it doesn't get
  # set below. (typically, it should be set from computing the walltime and
  # calling assign_job_queue() to get the corresponding queue type)
  destPartition <- "LocalQ"

  if (phylo == "FALSE") {
    # If not phylogenetic analysis, split up the sequences, run blast and full analysis
    num_seqs <- get_sequences(sequences, dir_path = dir, separate = TRUE)

    # calculate estimated walltime
    t_sec_estimate <- ifelse(
      type == "full" || type == "dblast",
      advanced_opts2est_walltime(advanced_options_names, num_seqs, n_hits = NHITS, verbose = TRUE),
      advanced_opts2est_walltime(advanced_options_names, num_seqs, verbose = TRUE)
    )

    # determine whether it's in the long or short queue
    # (we map the names 'long' and 'short' to our actual partition names,
    # "LocalQLong" and "LocalQ", respectively)
    destPartition <- ifelse(
      assign_job_queue(t_sec_estimate) == "long",
      "LocalQLong",
      "LocalQ"
    )
    cat(
      file=stderr(),
      stringr::str_glue(
        "submit_full(): estimated time {t_sec_estimate}, ",
        "queue {destPartition}\n\n"
      )
    )
    flush.console()

    destQoS <- ifelse(destPartition == "LocalQLong", "longjobs", "shortjobs")

    cmd_full <- stringr::str_glue(
      "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos={destQoS} --partition {destPartition} ",
      "--job-name {make_job_name(job_code, type)} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --array 1-{num_seqs} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb ",
      "input.txt {DB} {NHITS} {EVAL} F {type}"
    )
    job_id <- submit_and_log(cmd_full)
    job_ids <- c(job_ids, job_id)
    num_runs <- num_runs + num_seqs
  } else {
    get_sequences(sequences, dir_path = dir, separate = FALSE)
  }

  # run analysis on query regardless of selected advanced options
  destPartitionQuery <- "LocalQ" # query run always goes to short queues
  cmd_query <- stringr::str_glue(
    "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos=shortjobs --partition {destPartitionQuery} ",
    "--job-name {make_job_name(job_code, paste0(type, '_query'))} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --time=27:07:00 ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb ",
    "{sequences} {DB} {NHITS} {EVAL} T {type}"
  )
  job_id <- submit_and_log(cmd_query)
  job_ids <- c(job_ids, job_id)

  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")

  # schedule a 'summary' job to run when all the job_ids have completed or failed
  # (note that if no email was supplied, the job won't be scheduled at all)
  submit_summary_job(job_ids, submitter_email, dir, job_code)
}

submit_blast <- function(dir = "/data/scratch", blast = "~/test.fa", seqs = "~/seqs.fa", ncbi = FALSE, job_code=NULL, submitter_email=NULL, advanced_options=NULL, get_slurm_mails=FALSE) {
  # starts analysis for an input blast tsv
  # a query sequence(s) file can be provided,
  # or the sequences can be parsed from the Query column of input blast table
  setwd(dir)

  advanced_options_names <- names(advanced_options[advanced_options==TRUE])

  # write job submission params to file
  job_args <- list(
    submission_type = "blast",
    includes_ncbi_acc = ncbi,
    submitter_email = submitter_email,
    advanced_options = advanced_options_names,
    job_code = job_code
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  # keep track of the job ids so we can join on them for the summary job
  job_ids <- c()

  df_blast <- read_tsv(blast, col_names = web_blastp_hit_colnames)

  # get a df of only the unique query protein(s)
  df_query <- df_blast %>% distinct(Query, .keep_all = TRUE)
  # !!! 
  # the AccNum column here is overwritten for downstream joining between ipr-da columns and 
  # blast data ipr2da.R's append_ipr()
  # AccNum is overwritten by Query only for quick analyses on query proteins. 
  # this is to populate the Query data tab and not the Homolog data tab.
  # the intermediate query file currently, wrongly, carries the homology information as well
  # but this is fixed downstream in the app in molevolvr_app/scripts/MolEvolData_class.R 
  # where the spurious columns are removed before displaying in the web-app. 
  # These cols need to be removed upstream (before intermediate files are written) at a later time. 
  # #NextMilestone
  # !!!
  df_query$AccNum <- df_query$Query
  write_tsv(df_query, "blast_query.tsv", col_names = FALSE)
  # write the col of unique accession numbers of queries to a file
  write(df_query$Query, "accs.txt")

  # setup logfile table
  write(
    paste0(
      "START_DT\tSTOP_DT\tquery\tacc2info\tacc2fa\tcln_blast\tblast_clust\t",
      "clust2table\tiprscan\tipr2lineage\tipr2da\trps_blast\trps2da\tduration"
    ),
    file = "logfile.tsv"
  )

  # submit job for query proteins only
  destPartitionQuery <- "LocalQ" # query run always goes into short queue
  if (ncbi) {
    cmd_blast_query_ncbi <- stringr::str_glue(
      "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos=shortjobs --partition {destPartitionQuery} ",
      "--job-name {make_job_name(job_code, 'blast_query_ncbi')} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
      "blast_query.tsv T F"
    )
    job_id <- submit_and_log(cmd_blast_query_ncbi)
    job_ids <- c(job_ids, job_id)
  } else {
    cmd_blast_query <- stringr::str_glue(
      "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos=shortjobs --partition {destPartitionQuery} ",
      "--job-name {make_job_name(job_code, 'blast_query')} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
      "blast_query.tsv T T"
    )
    job_id <- submit_and_log(cmd_blast_query)
    job_ids <- c(job_ids, job_id)
  }

  # split dataframes by query and write them to their 
  # own output sub-directory
  for (df_split in split(df_blast, f = ~ Query)) {
    folder <-  paste0(unique(df_split$Query), "_blast")
    file <- paste0(unique(df_split$Query), ".dblast.tsv")
    dir.create(folder)
    write_tsv(df_split, file.path(folder, file), col_names = FALSE)
  }

  # calculate estimated walltime
  t_sec_estimate <- advanced_opts2est_walltime(
    advanced_options_names,
    n_inputs = nrow(df_blast),
    verbose = TRUE
  )
  # determine whether it's in the long or short queue
  # (we map the names 'long' and 'short' to our actual partition names,
  # "LocalQLong" and "LocalQ", respectively)
  destPartition <- ifelse(
    assign_job_queue(t_sec_estimate) == "long",
    "LocalQLong",
    "LocalQ"
  )
  cat(
    file=stderr(),
    stringr::str_glue(
      "submit_blast(): estimated time {t_sec_estimate}, ",
      "queue {destPartition}\n\n"
    )
  )
  flush.console()

  destQoS <- ifelse(destPartition == "LocalQLong", "longjobs", "shortjobs")

  # submit job array that is batched by the number of queries
  # for example,
  #   a blast table that only had one query protein will
  #   just be one job, but 2 queries protein would be split into 2 jobs, etc.
  # notably, the analysis on the query protein itself is still
  # done in the first submission above; a separate job
  n_queries <- nrow(df_query)
  cmd_blast <- stringr::str_glue(
    "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos={destQoS} --partition {destPartition} ",
    "--job-name {make_job_name(job_code, 'blast')} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --array 1-{n_queries} --time=27:07:00 ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
    "accs.txt F F"
  )
  job_id <- submit_and_log(cmd_blast)
  job_ids <- c(job_ids, job_id)

  # 1 is from the initial job on just the query proteins themselves
  # n_queries represents the total number of unique query proteins
  num_runs <- 1 + n_queries
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")

  # schedule a 'summary' job to run when all the job_ids have completed or failed
  # (note that if no email was supplied, the job won't be scheduled at all)
  submit_summary_job(job_ids, submitter_email, dir, job_code)
}

submit_ipr <- function(dir = "/data/scratch", ipr = "~/test.fa", seqs = "seqs.fa", ncbi = FALSE, blast = FALSE, DB = Sys.getenv("BLAST_DB", unset = "refseq"), NHITS = Sys.getenv("BLAST_HITS", unset = 100), EVAL = Sys.getenv("BLAST_EVALUE", unset = 0.00001), job_code=NULL, submitter_email=NULL, advanced_options=NULL, get_slurm_mails=FALSE) {
  setwd(dir)

  advanced_options_names <- names(advanced_options[advanced_options==TRUE])

  # write job submission params to file
  job_args <- list(
    submission_type = "interproscan",
    homology_search = blast,
    database = ifelse(blast == FALSE, NA, DB), # only include evalue, DB, & NHITS for blast jobs
    nhits = ifelse(blast == FALSE, NA, NHITS),
    includes_ncbi_acc = ncbi,
    submitter_email = submitter_email,
    advanced_options = advanced_options_names,
    job_code = job_code
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  job_ids <- c()

  # initialize counter for job progress
  num_runs <- 0
  # setup log file column headers
  write(
    paste0(
      "START_DT\tSTOP_DT\tquery\tdblast\tacc2info\t",
      "dblast_cleanup\tacc2fa\tblast_clust\tclust2table\t",
      "iprscan\tipr2lineage\tipr2da\tduration"
    ), 
    file = "logfile.tsv"
  )
  ipr_in <- read_tsv(ipr, col_names = TRUE)
  queries <- unique(ipr_in$AccNum)

  # calculate estimated walltime
  t_sec_estimate <- ifelse(
    blast,
    advanced_opts2est_walltime(advanced_options_names, length(queries), n_hits = NHITS, verbose = TRUE),
    advanced_opts2est_walltime(advanced_options_names, length(queries), verbose = TRUE)
  )
  # determine whether it's in the long or short queue
  # (we map the names 'long' and 'short' to our actual partition names,
  # "LocalQLong" and "LocalQ", respectively)
  destPartition <- ifelse(
    assign_job_queue(t_sec_estimate) == "long",
    "LocalQLong",
    "LocalQ"
  )
  cat(
    file=stderr(),
    stringr::str_glue(
      "submit_iprscan(): estimated time {t_sec_estimate}, ",
      "queue {destPartition}\n\n"
    )
  )
  flush.console()

  destQoS <- ifelse(destPartition == "LocalQLong", "longjobs", "shortjobs")

  if (ncbi) {
    # get seqs from accession number column
    acc2fa(queries, outpath = "seqs.fa")
  }
  if (blast) {
    # if blast separate the query sequences and do blast+full analysis
    seq_count <- get_sequences(seqs, dir_path = dir, separate = TRUE)
    num_runs <- num_runs + seq_count
    cmd_ipr_homology <- stringr::str_glue(
      "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos={destQoS} --partition {destPartition} ",
      "--job-name {make_job_name(job_code, 'ipr_homology')} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --array 1-{seq_count} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb ",
      "input.txt F T {DB} {NHITS} {EVAL}"
    )
    job_id <- submit_and_log(cmd_ipr_homology)
    job_ids <- c(job_ids, job_id)
  } else {
    write(queries, "accs.txt")
  }
  # add the query job to total runs
  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
  destPartitionQuery <- "LocalQ" # query run always goes to short queue
  # always do analysis on interpro file
  cmd_ipr_query <- stringr::str_glue(
    "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos=shortjobs --partition {destPartitionQuery} ",
    "--job-name {make_job_name(job_code, 'ipr_query')} --output=slurm_%x_%j.out --error=slurm_%x_%j.err --time=27:07:00 ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb ",
    "{ipr} T F"
  )
  job_id <- submit_and_log(cmd_ipr_query)
  job_ids <- c(job_ids, job_id)

  # schedule a 'summary' job to run when all the job_ids have completed or failed
  # (note that if no email was supplied, the job won't be scheduled at all)
  submit_summary_job(job_ids, submitter_email, dir, job_code)
}

submit_split_by_domain <- function(
  dir, sequences, DB = "refseq", NHITS = 5000,
  EVAL = 0.0001, phylo = "FALSE", type = "full",
  job_code = NULL, submitter_email = NULL, advanced_options = NULL,
  get_slurm_mails = FALSE
) {
  setwd(dir)
  # write job submission params to file
  job_args <- list(
    submission_type = type,
    database = ifelse(phylo == FALSE, DB, NA),
    nhits = ifelse(phylo == FALSE, NHITS, NA),
    evalue = ifelse(phylo == FALSE, EVAL, NA),
    submitter_email = submitter_email,
    advanced_options = advanced_options
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")
  # the pre-processing is not expensive 
  # (just an interproscan run on a FASTA and post-hoc data wrangling)
  # so we assing this job to the short queue always
  destQoS <- "shortjobs"
  destPartition <- "LocalQ"
  cmd_split_by_domain <- stringr::str_glue(
    "sbatch {make_email_args(submitter_email, get_slurm_mails)} --qos={destQoS} --partition ",
    "{destPartition} --job-name {make_job_name(job_code, 'fa2domain')} --time=27:07:00 ",
    "--output=split_by_domain-slurm-%j.out --error=split_by_domain-slurm-%j.err ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/split_by_domain-runner.R ",
    "{dir} {sequences} {DB} {NHITS} {EVAL} {phylo} {type} {job_code} {submitter_email} {advanced_options}"
  )
  submit_and_log(cmd_split_by_domain)
}
