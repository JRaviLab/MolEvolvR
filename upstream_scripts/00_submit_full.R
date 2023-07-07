suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
options(readr.show_col_types = FALSE) # silence read tsv col types
library(yaml)

get_sequences <- function(sequences, acc_file_path = "accs.txt", dir_path = "~", separate = TRUE) {
  seqs <- readAAStringSet(sequences)
  cln_names <- c()
  for (accnum in names(seqs)) {
    if (grepl("\\|", accnum)) {
      accnum_cln <- strsplit(accnum, "\\|")[[1]][2]
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

 #####     FA: i've left my debugging code here as an example of how to write messages
 #####     to the docker container log and actually be able to see them
 #####     (without flush.console(), you have to wait until the output buffer is full,
 #####     which in my recent experience didn't happen until i restarted the container)
submit_and_log <- function(cmd) {
  submit_result <- system(cmd)
  cat(file=stderr(), paste0("Ran command ", cmd, ", result: ", submit_result))
  flush.console()
}

submit_full <- function(dir = "/data/scratch", DB = "refseq", NHITS = 5000, EVAL = 0.0005, sequences = "~/test.fa", phylo = "FALSE", by_domain = "FALSE", domain_starting = "~/domain_seqs.fa", type = "full", job_code=NULL) {
  setwd(dir)
  # write job submission params to file
  job_args <- list(
    submission_type = type,
    database = ifelse(phylo == FALSE, DB, NA),
    nhits = ifelse(phylo == FALSE, NHITS, NA),
    evalue = ifelse(phylo == FALSE, EVAL, NA)
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  num_runs <- 0
  write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")
  if (phylo == "FALSE") {
    # If not phylogenetic analysis, split up the sequences, run blast and full analysis
    num_seqs <- get_sequences(sequences, dir_path = dir, separate = TRUE)
    cmd_full<- paste0("qsub -N ",
      make_job_name(job_code, type), " -t 1-", num_seqs,
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb 'input.txt ",
      DB, " ", NHITS, " ", EVAL, " F ", type, "'"
    )
    submit_and_log(cmd_full)
    num_runs <- num_runs + num_seqs
  } else {
    get_sequences(sequences, dir_path = dir, separate = FALSE)
  }
  # do analysis on query regardless of selected analysis
  if (by_domain == "TRUE") {
    submit_cmd_full_by_domain <- paste0("qsub -N ", 
      make_job_name(job_code, type),
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb '",
      domain_starting, " ", DB, " ", NHITS, " ", EVAL," T ", type, "'"
    )
    submit_and_log(submit_cmd_full_by_domain)
  } else {
    cmd_full <- paste0("qsub -N ", make_job_name(job_code, type), 
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb '", 
      sequences, " ", DB, " ", NHITS, " ", EVAL," T ", type, "'"
    )
    submit_and_log(cmd_full)
  }
  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
}

submit_blast <- function(dir = "/data/scratch", blast = "~/test.fa", seqs = "~/seqs.fa", ncbi = FALSE, job_code=NULL) {
  # starts jobs for BLAST output app submissions
  # a query sequence(s) file can be provided,
  # or the sequences can be parsed from the AccNum column of input blast table
  setwd(dir)

  # write job submission params to file
  job_args <- list(
    submission_type = "blast",
    includes_ncbi_acc = ncbi
  )
  yml <- yaml::as.yaml(job_args)
  write(yml, "job_args.yml")

  df_blast <- read_tsv(blast, col_names = web_blastp_hit_colnames)
  # get the accesion numbers (accnums) for homologs only
  homolog_accnums <- df_blast %>% filter(AccNum != Query) %>% select(AccNum) %>% distinct()
  n_homologs <- nrow(homolog_accnums)
  write_tsv(homolog_accnums, "accs.txt", col_names = FALSE)
  # make a folder for each homologous protein and write the the filtered rows
  for (homolog in base::unique(df_blast$AccNum)) {
    folder <- paste0(homolog, "_blast")
    file <- paste0(homolog, ".dblast.tsv")
    df_homolog <- df_blast %>% filter(AccNum == homolog)
    dir.create(folder)
    write_tsv(df_homolog, file.path(folder, file), col_names = FALSE)
  }
  cmd_blast_homologs <- paste0(
    "qsub -N ", make_job_name(job_code, "blast_homologs"), 
    " -t 1-", n_homologs, 
    " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb 'accs.txt", " F F'"
  )
  submit_and_log(cmd_blast_homologs)

  # select rows for blast query proteins (note: this is atypical, most submissions will only have a single query protein)
  # using 2 methods
  # 1. if the accnums match between query and hit, then use that
  df_query <- df_blast %>% 
    filter(Query == AccNum) %>% 
    distinct(Query, .keep_all = TRUE) # ensure no duplicates w/ distinct()
  # 2. if there's no match, then use the highest PcIdentity record 
    # first, find which queries did not have a hit w/ a matching accnum
  remaining_queries <- setdiff(df_blast$Query, df_query$Query)
    # second, sort by PcIdentity and use these records as the closest match 
    # to the query
  remaining_query_rows <- df_blast %>%
    filter(Query %in% remaining_queries) %>%
    arrange(desc(PcIdentity)) %>%
    distinct(Query, .keep_all = TRUE)
  # resulting table containing only rows w/ unique query accnums
  df_query <- bind_rows(df_query, remaining_query_rows)
  write_tsv(df_query, "blast_query.tsv", col_names = FALSE)

  # setup logfile table
  write("START_DT\tSTOP_DT\tquery\tacc2info\tacc2fa\tcln_blast\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\trps_blast\trps2da\tduration", "logfile.tsv")
  # do analysis on the queries
  if (ncbi) {
    cmd_blast_query_ncbi <- paste0("qsub -N ", make_job_name(job_code, "blast_query_ncbi"), 
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb 'blast_query.tsv ",
      " T F'"
    )
    submit_and_log(cmd_blast_query_ncbi)
  } else {
    cmd_blast_query <- paste0("qsub -N ", make_job_name(job_code, "blast_query"),
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb 'blast_query.tsv ",
      " T T'"  
    )
    submit_and_log(cmd_blast_query)
  }
  # n homologs + query
  num_runs <- 1 + n_homologs 
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
}

submit_ipr <- function(dir = "/data/scratch", ipr = "~/test.fa", seqs = "seqs.fa", ncbi = FALSE, blast = FALSE, DB = "refseq", NHITS = 5000, EVAL = 0.0005, job_code=NULL) {
  setwd(dir)

  # write job submission params to file
  job_args <- list(
    submission_type = "interproscan",
    homology_search = blast,
    database = ifelse(blast == FALSE, NA, DB), # only include evalue, DB, & NHITS for blast jobs
    nhits = ifelse(blast == FALSE, NA, NHITS),
    includes_ncbi_acc = ncbi
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
    cmd_ipr_homology <- paste0(
      "qsub -N ", make_job_name(job_code, "ipr_homology"), " -t 1-", seq_count, 
      " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb 'input.txt ", 
      "F T ", DB, " ", NHITS, " ", EVAL, "'"
    )
    submit_and_log(cmd_ipr_homology)
  } else {
    write(queries, "accs.txt")
  }
  # add the query job to total runs
  num_runs <- num_runs + 1
  write(paste0("0/", num_runs, " analyses completed"), "status.txt")
  # always do analysis on interpro file
  cmd_ipr_query <- paste0(
    "qsub -N ", make_job_name(job_code, "ipr_query"), 
    " /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb '", 
    ipr, " T F'"
  )
  submit_and_log(cmd_ipr_query)
}
