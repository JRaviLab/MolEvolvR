library("Biostrings")
args = commandArgs(trailingOnly = TRUE)
in_fa = readAAStringSet(args[1])
in_ipr = read.csv(args[2], sep="\t", header=FALSE)
out_file = args[3]
for (item in names(in_fa)){
  accession = unlist(strsplit(item, " "))[1]
  print(accession)
  for (i in 1:nrow(in_ipr)){
    domain_count = 0
    if (in_ipr[i, 1] == accession && (in_ipr[i,4] == "Pfam" || in_ipr[i,4] == "Gene3D")){
      sequence = toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0(">",accession, "_", in_ipr[i,4], "-", domain_count, " ", in_ipr[i,5], " ", in_ipr[i,6] )
      complete_seq = paste(header, sequence, sep="\n")
      cat(complete_seq, file=out_file, append=TRUE, sep = "\n")
      domain_count <- domain_count + 1
    }
    else{
    }
  }
}
