library("Biostrings")
args = commandArgs(trailingOnly = TRUE)
in_fa = readAAStringSet(args[1])
in_ipr = read.csv(args[2], sep="\t", header=FALSE)
out_file = args[3]
for (item in names(in_fa)){
  accession = unlist(strsplit(item, " "))[1]
  pfam_count = 0
  gene3d_count =0
  for (i in 1:nrow(in_ipr)){
    if (in_ipr[i, 1] == accession && (in_ipr[i,4] == "Pfam")){
      sequence = toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0(">",accession, "_", in_ipr[i,4], "-", pfam_count, " ", in_ipr[i,5], " ", in_ipr[i,6] )
      complete_seq = paste(header, sequence, sep="\n")
      cat(complete_seq, file=out_file, append=TRUE, sep = "\n")
      pfam_count <- pfam_count + 1
    }
    else if (in_ipr[i, 1] == accession && in_ipr[i,4] == "Gene3D"){
      sequence = toString(subseq(in_fa[item], in_ipr[i, 7], in_ipr[i, 8]))
      header <- paste0(">",accession, "_", in_ipr[i,4], "-", gene3d_count, " ", in_ipr[i,5], " ", in_ipr[i,6] )
      complete_seq = paste(header, sequence, sep="\n")
      cat(complete_seq, file=out_file, append=TRUE, sep = "\n")
      gene3d_count <- gene3d_count + 1
    }
    else{
    }
  }
}
