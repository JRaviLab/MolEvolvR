## ############ ##
## COLUMN names ##
## ############ ##
## BLAST
## Web-BLAST
web_blastp_colnames <- c("Query", "AccNum",
                         "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                         "QStart", "QEnd", "SStart", "SEnd",
                         "EValue", "BitScore", "PcPosOrig")

web_blastx_colnames <- c("Query", "AccNum",
                        "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                        "QStart", "QEnd", "SStart", "SEnd",
                        "EValue", "BitScore", "PcPosOrig",
                        "QSFrames") # specific to "blastx"
## Commandline BLAST
# cl_blast_colnames_orig <- c("qacc", "sacc", "sseqid",
#                             "sallseqid", "stitle", "sscinames", "staxids",
#                             "pident", "length", "mismatch", "gapopen",
#                             "qstart", "qend", "qlen",
#                             "sstart", "send", "slen",
#                             "evalue", "bitscore", "ppos")
# cl_blast_colnames_orig_renamed <- c("Query_AccNum", "sacc", "AccNum",
#                             "sallseqid", "stitle", "sscinames", "staxids",
#                             "pident", "length", "mismatch", "gapopen",
#                             "qstart", "qend", "qlen", "sstart", "send",
#                             "evalue", "bitscore", "ppos", "slen",
#                             "ppos_adjusted", "ClusterID")

cl_blast_colnames <- c("Query", "SAccNum", "AccNum",
                       "SAllSeqID", "STitle", "Species", "TaxID",
                       "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                       "QStart", "QEnd", "QLength",
                       "SStart", "SEnd", "SLength",
                       "EValue", "BitScore", "PcPosOrig",
                       "PcPositive", "ClusterID")  # post-cleanup

## IPRSCAN
# ipr_colnames_orig <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
#                        "DB_ID", "SignDesc", "StartLoc", "StopLoc", "Score",
#                        "Status", "RunDate", "IPRAcc", "IPRDesc")

ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                  "Status", "RunDate", "IPRAcc", "IPRDesc")

## RPSBLAST
# rps_colnames_orig <- c("qacc", "sacc", "sseqid", "pident", "ppos",
#                        "length", "mismatch",
#                        "qstart", "qend", "sstart", "send",
#                        "evalue", "bitscore", "staxids")
# rps_colnames_orig_renamed <- c("AccNum", "ID", "sseqid",
#                                "pident", "length", "mismatch",
#                                "qstart", "qend", "sstart", "send",
#                                "evalue", "bitscore")
rps_colnames <- c("AccNum", "DB.ID", "DBSeqID",
                  "PcIdentity.Dom", "PcPosOrig.Dom", #"PcPos.Dom", # Ppos missing
                  "AlnLength", "Mismatch",
                  "SStart", "SEnd", "DStart", "DEnd",
                  "EValue", "BitScore", "TaxID") # TaxID missing (NA); remove?

## IPG and Lineage maps
# ipg_colnames_orig <- c("Id", "Source", "Nucleotide Accession",
#                        "Start", "Stop", "Strand"  ,
#                        "Protein", "Protein Name",
#                        "Organism", "Strain", "Assembly")
ipg_colnames <- c("IPG.ID", "Source", "NucAccNum",
                  "NucStart", "NucStop", "Strand",
                  "AccNum", "ProtDesc",
                  "Species", "SppStrain", "AssemblyID")

## Final combined data frame that loads on the webapp
combo_colnames <- c("Query", "AccNum", "Species", "TaxID", "Lineage",
                    "PcPositive", "ClusterID",
                    # "Leaf", # MISSING (useful for all dataviz)
                    # "AssemblyID", "GeneName", "ProtDesc", # MISSING NOW!?!
                    "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
                    "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP",
                    "DomArch.SMART", "DomArch.TIGR")
