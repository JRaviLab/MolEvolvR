library("Biostrings")
in_fa = readAAStringSet("~/Downloads/example2.fa")
in_ipr = read.csv('~/Desktop/iprscan5-R20210528-185016-0482-14245416-p1m.tsv', sep="\t",
                  header=FALSE)
out_file = "~/Downloads/sep_fa.fa"
for (item in names(in_fa)){
  for (i in 1:nrow(in_ipr)){
    if (in_ipr[i, 1] == item){
      sequence = toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0("< ",item, " ", in_ipr[i,5], " ", in_ipr[i,6] )
      complete_seq = paste(header, sequence, sep="\n")
      cat(complete_seq, file=out_file, append=TRUE, sep = "\n")
    }
    else{
    }
  }
}

