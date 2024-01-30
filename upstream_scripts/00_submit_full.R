suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
options(readr.show_col_types = FALSE, readr.show_types = FALSE) # silence read tsv col types
library(yaml)

# source assign_job_queue.R via apps and scripts common ancestor
# assign_job_queue.R contains functions for computing how long jobs are estimated to take 
common_root <- rprojroot::has_file(".molevol_root")
source(common_root$find_file("molevol_scripts", "R", "assign_job_queue.R"))

get_sequences <- function(sequences, acc_file_path = "accs.txt", dir_path = "~", separate = TRUE) {
  seqs <- readAAStringSet(sequences)
  cln_names <- c()
  for (accnum in names(seqs)) {
    if (grepl("\\|", accnum)) {
      accnum_cln <- strsplit(accnum, "\\|")[[1]][2]
      accnum_cln <- stringr::str_split_1(accnum_cln, " ")[1]
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

# print info on slurm submission errors
submit_and_log <- function(cmd, exit = FALSE) {
  submit_result <- system(cmd)
  if (submit_result != 0L && exit == TRUE) stop(paste0("failed to submit job; error code: ", submit_result, "\n", "cmd=", cmd))
  cat(file=stderr(), paste0("Ran command ", cmd, ", result: ", submit_result, "\n"))
  flush.console()
}

submit_full <- function(dir = "/data/scratch", DB = "refseq", NHITS = 5000, EVAL = 0.0005, sequences = "~/test.fa", phylo = "FALSE", by_domain = "FALSE", domain_starting = "~/domain_seqs.fa", type = "full", job_code=NULL, submitter_email=NULL, advanced_options=NULL) {
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
    advanced_options = advanced_options_names
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  num_runs <- 0
  write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")

  # assign destPartition to a known partition name just in case it doesn't get
  # set below. (typically, it should be set from computing the walltime and
  # calling assign_job_queue() to get the corresponding queue type)
  destPartition <- "LocalQ"

  if (phylo == "FALSE") {
    # If not phylogenetic analysis, split up the sequences, run blast and full analysis
    num_seqs <- get_sequences(sequences, dir_path = dir, separate = TRUE)

    # compute estimated walltime, so we can store that in the collection
    t_sec_esimate <- input_opts2est_walltime(advanced_options_names, num_seqs)

    # determine whether it's in the long or short queue
    # (we map the names 'long' and 'short' to our actual partition names,
    # "LocalQLong" and "LocalQ", respectively)
    destPartition <- ifelse(assign_job_queue(t_sec_esimate) == "long", "LocalQLong", "LocalQ")
    cat(file=stderr(), stringr::str_glue("submit_full(): estimated time {t_sec_esimate}, queue {destQueue}\n\n"))
    flush.console()

    cmd_full <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, type)} --array 1-{num_seqs} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb ",
      "input.txt {DB} {NHITS} {EVAL} F {type}"
    )
    submit_and_log(cmd_full)
    num_runs <- num_runs + num_seqs
  } else {
    get_sequences(sequences, dir_path = dir, separate = FALSE)
  }

  # do analysis on query regardless of selected analysis
  if (by_domain == "TRUE") {
    cmd_by_domain <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, type)} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb ",
      "{domain_starting} {DB} {NHITS} {EVAL} T {type}"
    )
    submit_and_log(cmd_by_domain)
  } else {
    cmd_query <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, type)} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb ",
      "{sequences} {DB} {NHITS} {EVAL} T {type}"
    )
    submit_and_log(cmd_query)
  }

  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
}

submit_blast <- function(dir = "/data/scratch", blast = "~/test.fa", seqs = "~/seqs.fa", ncbi = FALSE, job_code=NULL, submitter_email=NULL, advanced_options=NULL) {
  # starts analysis for an input blast tsv
  # a query sequence(s) file can be provided,
  # or the sequences can be parsed from the Query column of input blast table
  setwd(dir)

  # write job submission params to file
  job_args <- list(
    submission_type = "blast",
    includes_ncbi_acc = ncbi,
    submitter_email = submitter_email
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

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

  # as in submit_full(), we assign a default partition in case one wasn't computed
  destPartition <- "LocalQ"

  # submit job for query proteins only
  if (ncbi) {
    cmd_blast_query_ncbi <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, 'blast_query_ncbi')} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
      "blast_query.tsv T F"
    )
    submit_and_log(cmd_blast_query_ncbi)
  } else {
    cmd_blast_query <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, 'blast_query')} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
      "blast_query.tsv T T"
    )
    submit_and_log(cmd_blast_query)
  }

  # split dataframes by query and write them to their 
  # own output sub-directory
  for (df_split in split(df_blast, f = ~ Query)) {
    folder <-  paste0(unique(df_split$Query), "_blast")
    file <- paste0(unique(df_split$Query), ".dblast.tsv")
    dir.create(folder)
    write_tsv(df_split, file.path(folder, file), col_names = FALSE)
  }

  # submit job array that is batched by the number of queries
  # for example, 
  #   a blast table that only had one query protein will
  #   just be one job, but 2 queries protein would be split into 2 jobs, etc.
  # notably, the analysis on the query protein itself is still 
  # done in the first submission above; a separate job
  n_queries <- df_blast %>% select(Query) %>% distinct() %>% nrow()
  cmd_blast_homologs <- stringr::str_glue(
    "sbatch --partition {destPartition} --job-name {make_job_name(job_code, 'blast')} --array 1-{n_queries} --time=27:07:00 ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb ",
    "accs.txt F F"
  )
  submit_and_log(cmd_blast_homologs)

  # 1 is from the initial job on just the query proteins themselves
  # n_queries represents the total number of unique query proteins
  num_runs <- 1 + n_queries
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
}

submit_ipr <- function(dir = "/data/scratch", ipr = "~/test.fa", seqs = "seqs.fa", ncbi = FALSE, blast = FALSE, DB = "refseq", NHITS = 5000, EVAL = 0.0005, job_code=NULL, submitter_email=NULL, advanced_options=NULL) {
  setwd(dir)

  # write job submission params to file
  job_args <- list(
    submission_type = "interproscan",
    homology_search = blast,
    database = ifelse(blast == FALSE, NA, DB), # only include evalue, DB, & NHITS for blast jobs
    nhits = ifelse(blast == FALSE, NA, NHITS),
    includes_ncbi_acc = ncbi,
    submitter_email = submitter_email
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

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

  # as in submit_full(), we assign a default partition in case one wasn't computed
  destPartition <- "LocalQ"

  ipr_in <- read_tsv(ipr, col_names = TRUE)
  queries <- unique(ipr_in$AccNum)
  if (ncbi) {
    # get seqs from accession number column
    acc2fa(queries, outpath = "seqs.fa")
  }
  if (blast) {
    # if blast separate the query sequences and do blast+full analysis
    seq_count <- get_sequences(seqs, dir_path = dir, separate = TRUE)
    num_runs <- num_runs + seq_count
    cmd_ipr_homology <- stringr::str_glue(
      "sbatch --partition {destPartition} --job-name {make_job_name(job_code, 'ipr_homology')} --array 1-{seq_count} --time=27:07:00 ",
      "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb ",
      "input.txt F T {DB} {NHITS} {EVAL}"
    )
    submit_and_log(cmd_ipr_homology)
  } else {
    write(queries, "accs.txt")
  }
  # add the query job to total runs
  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
  # always do analysis on interpro file
  cmd_ipr_query <- stringr::str_glue(
    "sbatch --partition {destPartition} --job-name {make_job_name(job_code, 'ipr_query')} --time=27:07:00 ",
    "/data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb ",
    "{ipr} T F"
  )
  submit_and_log(cmd_ipr_query)
}
