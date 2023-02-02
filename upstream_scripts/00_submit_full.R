library(tidyverse)
library(Biostrings)

get_sequences <- function(sequences, acc_file_path = "accs.txt", dir_path = "~", separate = TRUE){

    test_valid_accession <- function(header) {
        # checks FASTA input for valid accession
        # accession number prefix format guide
        # https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
        # https://academic.oup.com/nar/article/33/suppl_1/D501/2505241?login=false
        # patterns: 
        # 1 GenBank
        # 2 Swiss-Prot/UniProtKB
        # 3 RefSeq
        pattern <- paste("[A-Z]{3}[0-9]{5}|[A-Z]{3}[0-9]{7}",
                         "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
                         "[A-Z][A-Z]_[0-9]{6}|[A-Z][A-Z]_[0-9]{8}|[A-Z][A-Z]_[0-9]{9}|[A-Z][A-Z]_[A-Za-z]{4}",
                         sep="|")
        result <- grepl(pattern, header)
        return(result)
    }

    seqs <- readAAStringSet(sequences)
    for (seq_num in 1:length(seqs)) {
        header <- names(seqs[seq_num])
        # clean header
        header <- gsub("[^A-Za-z]", "", header)
        if (test_valid_accession(header) == FALSE) {
            system('mkdir seq2acc_data')
            # use tempfile to prevent duplicate headers overwriting the same file
            single_fa_path <- tempfile(header,
                                       tmpdir=paste0(getwd(), '/', 'fa2acc_data'),
                                       fileext = '.fa')
            writeXStringSet(seqs[seq_num], single_fa_path)
            # run BLAST
            system(paste('sh',
                         '/data/research/jravilab/molevol_scripts/upstream_scripts/seq2acc.sh',
                         single_fa_path)
            )
             # accessing blast output
                # blast output is saved at PATH=single_fa_path
            blast_output_file <- paste0(strsplit(single_fa_path, '.fa'), '.csv')
            df <- read.csv(blast_output_file, header=FALSE)
            # if an entry has 100% seq similarity (col3), get the accnum in col2
            accnum <- df[df[,3] == 100][2]
            # change the header to the accnum found from blast
            names(seqs)[seq_num] <- accnum
        }
    }

    cln_names <- c()
    for (accnum in names(seqs)){
        if (grepl("\\|", accnum)){
            accnum_cln <- strsplit(accnum, "\\|")[[1]][2]
        }
        else{
            accnum_cln <- strsplit(accnum, " ")[[1]][1]
        }
        cln_names <- append(cln_names, accnum_cln)
        write(accnum_cln, file = acc_file_path, append = TRUE)
        if (separate){
        write(paste0(dir_path, "/", accnum_cln, ".faa"), file = "input.txt", append = TRUE)
        write(paste0(">", accnum_cln), file = paste0(accnum_cln, ".faa"), append = TRUE)
        write(toString(seqs[accnum]), file = paste0(accnum_cln, ".faa"), append = TRUE)
        }
    }
    names(seqs) <- cln_names
    writeXStringSet(seqs, sequences, format="fasta")
    return(length(seqs))
}

submit_full <- function(dir = "/data/scratch", DB = "refseq", NHITS = 5000, EVAL= 0.0005, sequences = "~/test.fa", phylo = "FALSE", by_domain = "FALSE", domain_starting = "~/domain_seqs.fa", type = "full"){
    setwd(dir)
    num_runs <- 0
    write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")
    if (phylo == "FALSE"){
        # If not phylogenetic analysis, split up the sequences, run blast and full analysis
        num_seqs <- get_sequences(sequences, dir_path = dir, separate = TRUE)
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F 'input.txt ", DB, " ", NHITS, " ", EVAL, " F ", type, "'", " -t 1-", num_seqs))
        num_runs <- num_runs + num_seqs
    }
    else{
        get_sequences(sequences, dir_path = dir, separate = FALSE)
    }
    # do analysis on query regardless of selected analysis
    if (by_domain == "TRUE"){
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F '", domain_starting, " ", DB, " ", NHITS, " ", EVAL," T ", type, "'"))
    }
    else{
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F '", sequences, " ", DB, " ", NHITS, " ", EVAL," T ",type, "'"))
    }
    num_runs <- num_runs + 1
    write(paste0("0/", num_runs, " analyses completed"), "status.txt")
}

submit_blast <- function(dir = "/data/scratch", blast = "~/test.fa", seqs = "~/seqs.fa", ncbi = FALSE){ 
    #Accepts sequence file, if none given then download the ncbi sequences
    setwd(dir)
    num_runs <- 0
    df <- read_tsv(blast, col_names = web_blastp_hit_colnames)
    df_query <- df %>% distinct(Query, .keep_all = TRUE)
    df_query$AccNum <- df_query$Query
    write_tsv(df_query, "blast_query.tsv", col_names = FALSE)
    write(df_query$AccNum, "accs.txt")
    write("START_DT\tSTOP_DT\tquery\tacc2info\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")
    # do analysis on the queries
    if (ncbi){
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'blast_query.tsv ", " T F'"))
    }
    else{
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'blast_query.tsv ", " T T'"))
    }
    num_runs <- num_runs + 1
    for (query in unique(df$Query)){
        system(paste0("mkdir ", query, "_blast"))
        df_filter <- df %>% filter(query == Query)
        write_tsv(df_filter, paste0(query, "_blast/", query, ".dblast.tsv"), col_names = FALSE)
    }
    num_runs <- num_runs + length(unique(df$Query))
    write(paste0("0/", num_runs, " analyses completed"), "status.txt")
    # do analysis on the homologs
     system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'accs.txt ", " F F'", " -t 1-", length(unique(df$Query))))
}

submit_ipr <- function(dir = "/data/scratch", ipr = "~/test.fa", seqs = "seqs.fa", ncbi = FALSE, blast = FALSE, DB = "refseq", NHITS = 5000, EVAL= 0.0005){
    setwd(dir)
    num_runs <- 0
    write("START_DT\tSTOP_DT\tquery\tacc2info\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration", "logfile.tsv")
    ipr_in <- read_tsv(ipr, col_names = TRUE)
    queries <- unique(ipr_in$AccNum)
    if (ncbi){
        acc2fa(queries, outpath = "seqs.fa")
    }
    if (blast){
        # if blast separate the query sequences and do blast+full analysis
        seq_count <- get_sequences(seqs, dir_path = dir, separate = TRUE)
        num_runs <- num_runs + seq_count
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb -F 'accs.txt ", "F T ", DB , " ", NHITS, " ", EVAL, "'", " -t 1-", seq_count))
    }
    else{
        write(queries, "accs.txt")
    }
    num_runs <- num_runs + 1
    write(paste0("0/", num_runs, " analyses completed"), "status.txt")
    # always do analysis on interpro file
    system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_ipr.sb -F '", ipr, " T F'"))
    
}

#args <- commandArgs(trailingOnly = TRUE)

## call function
#submit_full(args[1], args[2], args[3], args[4], args[5], args[6])
