####### Accession 2 Lineage Example + Setup ########
# This script contains all the steps to generate the necessary files to run
# 'acc2lin()'
# It also includes an example run

### Imports:
source("R/acc2lin.R")
source("R/GCA2Lins.R")
source("R/create_lineage_lookup.R")


##### 1. Create Lineage Lookup table #####
# The first step is to create the lineage lookup table that maps taxIDs to
# lineage
# The 'create_lineage_lookup()' function achieves this by cleaning up the
# rankedlineage.dmp file, which can be downloaded at
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# Path to the ranked lineage dump file
lineage_dmp = "data/rankedlineage.dmp"

# Path to store the resulting lineage lookup table
lineagelookup_path = "data/lineagelookup.txt"

create_lineage_lookup(lineage_file = lineage_dmp, outfile = lineagelookup_path,
                      taxonomic_rank = "phylum")


##### 2. Download the Assembly Accession file ####
# The next step is to download the assembly summary file which maps GCA to the
# tax id.
# The file also contains columnts for the species taxid and organism name

# Path to store the downloaded assembly summary
Assembly_path = paste0("data/acc_files/assembly_summary", Sys.Date(), ".txt")
DownloadAssemblySummary(outpath = Assembly_path)

##### 3. Running acc2lin #####
# After the previous 2 steps have been completed, the accession numbers
# can be mapped to GCA_ID, TaxID, Lineage, and a few more columns

accessions = read_tsv("data/rawdata_tsv/all_clean.txt" , col_names = T) %>%
  pull(AccNum)

# Pull 5,000 sample accession numbers
accessions = sample(accessions, size = 5000, replace = F)

# Path to write the results of the efetch run against the ipg database to
# NULL will not save the file anywhere
ipgout_path = NULL

s_time = Sys.time()
lins = acc2lin(accessions = accessions, assembly_path = Assembly_path,
        lineagelookup_path = lineagelookup_path, ipgout_path = ipgout_path )
e_time = Sys.time()
elapsed = e_time - s_time
print(elapsed)

view(lins)
