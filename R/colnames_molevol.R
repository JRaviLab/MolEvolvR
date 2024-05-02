## ############ ##
## COLUMN names ##
## ############ ##
## BLAST
############
## Web-BLAST
############
## Downloaded as HIT-TABLE csv
# BLASTP and related protein BLASTs
# web_blastp_hit_colnames <- c(
#   "Query", "AccNum",
#   "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
#   "QStart", "QEnd", "SStart", "SEnd",
#   "EValue", "BitScore", "PcPosOrig"
# )
# BLASTX
# web_blastx_colnames <- c(
#   "Query", "AccNum",
#   "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
#   "QStart", "QEnd", "SStart", "SEnd",
#   "EValue", "BitScore", "PcPosOrig",
#   "QSFrames"
# ) # specific to "blastx"

# ## Downloaded as Descriptions csv
# # BLASTP and related protein BLASTs
# web_blastp_desc_colnames <- c(
#   "Description", "Species", "CommonName", "TaxID",
#   "BitScore", "TotalScore",
#   "PcQCover", "EValue", "PcIdentity",
#   "SLen", "AccNum"
# )
# # Ref: https://ncbiinsights.ncbi.nlm.nih.gov/2020/11/23/blast-new-columns/
# # Description,	Scientific Name,	Common Name,	Taxid,
# # Max Score,	Total Score,
# # Query Cover,	E value,	Per. ident,
# # Acc. Len	Accession


#####################
## Command line BLAST
#####################

# # pre-cleanup
# cl_blast_colnames <- c(
#   "Query", "SAccNum", "AccNum",
#   "SAllSeqID", "STitle", "Species", "TaxID",
#   "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
#   "QStart", "QEnd", "QLength",
#   "SStart", "SEnd", "SLength",
#   "EValue", "BitScore", "PcPosOrig"
# )

# post-cleanup
# cl_blast_postcln_cols <- c(
#   "Query", "AccNum",
#   "STitle", "Species", "TaxID", "Lineage", "Lineage_long", "Lineage_long_na", "Lineage_med", "Lineage_short",
#   "PcPositive", "PcIdentity", "AlnLength",
#   "SAccNum", "SAllSeqID",
#   "Mismatch", "GapOpen",
#   "QStart", "QEnd", "QLength",
#   "SStart", "SEnd", "SLength",
#   "EValue", "BitScore", "PcPosOrig", "QueryName"
# )
# cl_blast_postcln_cols <- c("Query", "SAccNum", "AccNum",
#                       "SAllSeqID", "STitle", "Species", "TaxID",
#                       "PcIdentity", "AlnLength", "Mismatch", "GapOpen",
#                       "QStart", "QEnd", "QLength",
#                       "SStart", "SEnd", "SLength",
#                       "EValue", "BitScore", "PcPosOrig",
# 			"PcPositive", "Lineage")


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


##########
## IPRSCAN
##########
# ipr_colnames_orig <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
#                        "DB_ID", "SignDesc", "StartLoc", "StopLoc", "Score",
#                        "Status", "RunDate", "IPRAcc", "IPRDesc")

# ipr_colnames <- c(
#   "AccNum", "SeqMD5Digest", "SLength", "Analysis",
#   "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
#   "Status", "RunDate", "IPRAcc", "IPRDesc"
# )

# post cleanup
##########################
## NEED TO BE REORDERED ##
##########################
# ipr_cln_colnames <- c(
#   "DB.ID", "TaxID", "AccNum.noV", "AccNum",
#   "SeqMD5Digest", "SLength", "Analysis", "SignDesc",
#   "StartLoc", "StopLoc", "Score", "Status", "RunDate",
#   "IPRAcc", "IPRDesc", "FullAccNum", "ProteinName",
#   "Length", "SourceDB", "Completeness", "Lineage",
#   "Species", "Name", "ShortName", "LookupTblDesc",
#   "ID", "Label"
# )

###########
## RPSBLAST
###########
# rps_colnames_orig <- c("qacc", "sacc", "sseqid", "pident", "ppos",
#                        "length", "mismatch",
#                        "qstart", "qend", "sstart", "send",
#                        "evalue", "bitscore", "staxids")
# rps_colnames_orig_renamed <- c("AccNum", "ID", "sseqid",
#                                "pident", "length", "mismatch",
#                                "qstart", "qend", "sstart", "send",
#                                "evalue", "bitscore")
# rps_colnames <- c(
#   "AccNum", "DB.ID", "DBSeqID",
#   "PcIdentity.Dom", "PcPosOrig.Dom", # "PcPos.Dom", # Ppos missing
#   "AlnLength", "Mismatch",
#   "SStart", "SEnd", "DStart", "DEnd",
#   "EValue", "BitScore", "TaxID"
# ) # TaxID missing (NA); remove?

#######################
## IPG and Lineage maps
#######################
# ipg_colnames_orig <- c("Id", "Source", "Nucleotide Accession",
#                        "Start", "Stop", "Strand"  ,
#                        "Protein", "Protein Name",
#                        "Organism", "Strain", "Assembly")
# ipg_colnames <- c(
#   "IPG.ID", "Source", "NucAccNum",
#   "NucStart", "NucStop", "Strand",
#   "AccNum", "Description",
#   "Species", "Spp.Strain", "AssemblyID"
# )

##################
## Assembly files
## Genbank, Refseq
##################
# assembly_colnames <- c(
#   "AssemblyID",
#   "bioproject", "biosample", "wgs_master", # not used
#   "RefseqCategory", "TaxID", "Spp.TaxID",
#   "Species", "Spp.Strain",
#   "isolate", "version_status", # not used
#   "assembly_level", "release_type", # not used
#   "GenomeStatus",
#   "seq_rel_date", "asm_name", "submitter", # not used
#   "AssemblyID.GBRS",
#   "paired_asm_comp", "ftp_path", # not used
#   "excluded_from_refseq", "relation_to_type_material"
# ) # not used
# assembly_sub_colnames <- c(
#   "TaxID", "Spp.TaxID", "Species", "Spp.Strain",
#   "RefseqCategory", "GenomeStatus",
#   "AssemblyID", "AssemblyID.GBRS"
# )

# assembly_colnames_orig <- c("assembly_accession", "bioproject", "biosample",
#  "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name",
#  "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type",
#  "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm",
#  "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material")

#################
## Lookup tables
## in common_data
#################
# lineage_lookup_colnames <- c("TaxID", "Species", "Lineage_long", "Lineage_long_na", "Lineage_med", "Lineage_short", "Lineage")
# domarch_lookup_colnames <- c("DB.ID", "ShortName", "Description", "ID")
# !! SC and LS will fix other piecemeal files based on these

######################
## FINAL UPLOADED DATA
######################
## Combined data frame that is loaded on to the webapp
# combo_colnames <- c(
#   "Query", "UID", "AccNum", "Species", "TaxID", "Lineage",
#   "PcPositive", "ClusterID", "QueryName",
#   # "AssemblyID", "GeneName", "Description", # MISSING NOW!?!
#   "DomArch.Pfam", "DomArch.COG", "DomArch.Gene3D",
#   "DomArch.TMHMM", "DomArch.Phobius", "DomArch.SignalP",
#   "DomArch.SMART", "DomArch.TIGR"
# )


################
## read tsv colnames
################
# lookup_table_cols <- cols(
#   DB.ID = col_character(),
#   ShortName = col_character(),
#   Description = col_character(),
#   ID = col_character()
# )

# iprscan_cols <- cols(
#   .default = col_character(),
#   TaxID = col_double(),
#   SLength = col_double(),
#   SignDesc = col_character(),
#   StartLoc = col_double(),
#   StopLoc = col_double(),
#   Score = col_double(),
#   Status = col_logical(),
#   IPRAcc = col_character(),
#   IPRDesc = col_character(),
#   Length = col_double(),
#   ShortName = col_character(),
#   LookupTblDesc = col_character(),
#   ID = col_character(),
#   Label = col_character()
# )

# ipr_cln_cols <- cols(
#   .default = col_character(),
#   TaxID = col_double(),
#   SLength = col_double(),
#   StartLoc = col_double(),
#   StopLoc = col_double(),
#   Score = col_double(),
#   Status = col_logical(),
#   IPRAcc = col_logical(),
#   IPRDesc = col_logical(),
#   Length = col_double(),
#   ID = col_logical()
# )

# lineage_map_cols <- c(
#   "double",
#   "character",
#   "character", "character", "character", "character", "character"
# )
