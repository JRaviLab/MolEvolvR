#library(tidyverse)
#library(Biostrings)

separate_sequences <- function(sequences, acc_file_path = "accs.txt", dir_path = "~"){
    seqs <- readAAStringSet(sequences)
    for (accnum in names(seqs)){
        if ("\\|" %in% accnum){
            accnum_cln <- strsplit(accnum, "\\|")[[1]][3]
        }
        else{
            accnum_cln <- strsplit(accnum, " ")[[1]][1]
        }
        write(accnum_cln, file = acc_file_path, append = TRUE)
        write(paste0(dir_path, "/", accnum_cln, ".faa"), file = "input.txt", append = TRUE)
        write(paste0(">", accnum_cln), file = paste0(accnum_cln, ".faa"), append = TRUE)
        write(toString(seqs[accnum]), file = paste0(accnum_cln, ".faa"), append = TRUE)
    }
    write(paste0("0/", (length(seqs) + 1), " analyses completed"), "status.txt")
    return(length(seqs))
}

submit_full <- function(dir = "/data/scratch", DB = "refseq", NHITS = 5000, EVAL= 0.0005, sequences = "~/test.fa", phylo = FALSE){
    setwd(dir)
    write("START_DT\tSTOP_DT\tquery\tdblast\tacc2info\tdblast_cleanup\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration\n", "logfile.tsv")
    if (!phylo){
        num_seqs <- separate_sequences(sequences, dir_path = dir)
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F 'input.txt ", DB, " ", NHITS, " ", EVAL, " F", "' -t 1-", num_seqs))
    }
    system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_full.sb -F '", sequences, " ", DB, " ", NHITS, " ", EVAL," T'"))
}

submit_blast <- function(dir = "/data/scratch", blast = "~/test.fa", seqs = "~/seqs.fa", ncbi = FALSE){ 
    #Accepts sequence file, if none given then download the ncbi sequences
    setwd(dir)
    df <- read_tsv(blast, col_names = web_blastp_hit_colnames)
    df_query <- df %>% distinct(Query, .keep_all = TRUE)
    df_query$AccNum <- df_query$Query
    write_tsv(df_query, "blast_query.tsv", col_names = FALSE)
         write("START_DT\tSTOP_DT\tquery\tacc2info\tacc2fa\tblast_clust\tclust2table\tiprscan\tipr2lineage\tipr2da\tduration\n", "logfile.tsv")
    if (ncbi){
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'blast_query.tsv ", " T F'"))
    }
    else{
        system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'blast_query.tsv ", " T T'"))
    }
    for (query in unique(df$Query)){
        system(paste0("mkdir ", query, "_blast"))
        df_filter <- df %>% filter(query == Query)
        write_tsv(df_filter, paste0(query, "_blast/", query, ".dblast.tsv"), col_names = FALSE)
    }
     system(paste0("qsub /data/research/jravilab/molevol_scripts/upstream_scripts/00_wrapper_blast.sb -F 'accs.txt ", " F F'", " -t 1-", length(unique(df$Query))))
}

#args <- commandArgs(trailingOnly = TRUE)

## call function
#submit_full(args[1], args[2], args[3], args[4], args[5], args[6])