setClass("MolEvolData",
                          slots = list(df = "data.frame",
                                        fasta_seq = "character",
                                        msa = "character",
                                        fasta_path = "character",
                                        msa_path = "character",
                                        ipr_path = "character",
                                        queries = "character"
                          )
                          )


process_wrapper_dir <- function(path)
{
  files = list.files(path = path)

  fasta_file = files[grepl("all_accnums\\.fa$",files)]
  ipr_file = files[grepl("iprscan\\.tsv", files)]
  cln_blast_file = files[grepl("refseq\\..+\\.cln.txt$", files)]

  fasta_path = paste0(path,"/", fasta_file)
  fasta_seq = read_file(fasta_path)

  ipr_path = paste0(path, "/", ipr_file)


  cln_blast = read_tsv(paste0(path,"/", cln_blast_file), col_names = T)
  # length(which(is.na(di$DomArch.SMART) == 1))

  wrapper_data <- new("MolEvolData", fasta_path = fasta_path, fasta_seq = fasta_seq,
                      ipr_path = ipr_path, df = cln_blast)

  return(wrapper_data)

}
