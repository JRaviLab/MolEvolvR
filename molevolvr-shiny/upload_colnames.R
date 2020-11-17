
# Web BLAST output
web_blast_colnames <- c("Query", "AccNum",
                        "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                        "QStart", "QEnd", "SStart", "SEnd",
                        "EValue", "BitScore", "PcPosOrig",
                        "QSFrames") # specific to "blastx"


# BLAST Command line
cl_blast_colnames <- c("Query", "SAccNum", "AccNum",
                       "SAllSeqID", "STitle", "Species", "TaxID",
                       "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
                       "QStart", "QEnd", "QLength",
                       "SStart", "SEnd", "SLength",
                       "EValue", "BitScore", "PcPosOrig",
                       "PcPositive", "ClusterID") # post-cleanup

# IPRSCAN (web+command-line)
ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                  "Status", "RunDate", "IPRAcc", "IPRDesc")

# RPSBLAST
rps_colnames <- c("AccNum", "DBID", "DBSeqID",
                  "PcIdentity", "PcPosOrig", # Ppos missing
                  "AlnLength", "Mismatch",
                  # Q here is Subject; S here is the matching domain; rename!
                  "SStart", "SEnd", "DStart", "DEnd",
                  "EValue", "BitScore", "TaxID") # TaxID missing (NA); remove?

# IPG
ipg_colnames <- c("IPG.ID", "Source", "NucAccNum",
                  "NucStart", "NucStop", "Strand",
                  "AccNum", "ProtDesc",
                  "Species", "SppStrain", "AssemblyID")
# Final ColNames
combo_colnames <- c("Query", "AccNum", "Species", "TaxID", "Lineage",
                    "PcPositive", "ClusterID",
                    # "Leaf", # MISSING (useful for all dataviz)
                    # "AssemblyID", "GeneName", "ProtDesc", # MISSING NOW!?!
                    "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
                    "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP")
