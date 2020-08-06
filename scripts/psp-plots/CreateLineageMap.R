library(tidyverse)
library(biomartr)

# ##### Use these to generate assembly summary files ########
# ## Over 680,000
# assembly_k_genbank <- getKingdomAssemblySummary("genbank")
#
# ## Over 200,000
# assembly_k_refseq <- getKingdomAssemblySummary("refseq")
#
# keep = c("assembly_accession", "taxid", "species_taxid", "organism_name")
# assembly_all <- bind_rows(assembly_k_genbank,assembly_k_refseq) %>%
#                   select(all_of(keep))
#
# write_tsv(assembly_all, "data/acc_files/assembly_summary20200706.txt")

#Goal:
lin_map <- read_tsv("data/acc_files/organisms2lineages_map_bae_20170828.txt")

# Go from the GCA_ID column to tax IDs using the assembly file
prot_data <- read_tsv("data/rawdata_tsv/all_clean.txt")
prot_data <- prot_data %>% select(AccNum, assembly_accession = GCA_ID)

mergedTax <- prot_data %>% left_join(assembly_all,by = "assembly_accession")
na_taxid <- which(is.na(mergedTax$taxid))

na_obs <- mergedTax[na_taxid,]

shorten_NA <- function(Lineage)
{
  first_NA = str_locate(Lineage, "NA")[1]
  if(is.na(first_NA) )
  {
    # No NAs
    # print(Lineage)
    return(Lineage)
  }
  else
  {
    shortened = substr(Lineage,1,(first_NA-1))
    return(shortened)
  }
}




# assembly_file <- read_tsv("C:/Users/samue/Google_Drive/GitHub/molevol/data/assembly_summary_genbank.20190312.txt", skip = 1
#                           , na = "na")
# assembly_file <- assembly_file %>% select(GCA_ID = "# assembly_accession", taxid, species_taxid, organism_name)
#
dont_include = c("X2", "X4", "X6", "X8","X10","X12", "X14","X16","X17", "X18","X20")
rankedLins <- read_tsv("C:/Users/samue/Google_Drive/GitHub/molevol/data/rankedlineage.dmp", col_names = F) %>%
  select(-all_of(dont_include)) %>% rename(taxid = X1)


merged_lins <- left_join(mergedTax, rankedLins, by = "taxid")

merged_lins <- merged_lins %>% unite(col = 'Lineage',X19:X5, sep = ">") %>%
                mutate(Lineage = map(Lineage,shorten_NA)) %>%
                unite(Lineage, Lineage, X3 )



# group <- rep(1:cl, length.out = nrow(taxID2lins_collapsed))
#
# taxID2lins_collapsed <- bind_cols(data.frame(group), taxID2lins_collapsed)
#
# cluster=new_cluster(cl); cluster_copy(cluster,c("map","shorten_NA", "str_locate","is.na","substr"))
#
# partitioned_group <- taxID2lins_collapsed %>% group_by(group) %>% partition(cluster)
#
#
# ### 2.334489 minutes Runtime: 4 cores
# start_time <- Sys.time()
# taxID2lins_noNA <- partitioned_group %>% mutate(Lineage2 = map(Lineage,shorten_NA)) %>% collect()
# end_time <- Sys.time()
# end_time - start_time





### 1.2175 minutes Runtime
# start_time2 <- Sys.time()
# taxID2lins_noNA <- taxID2lins_collapsed %>% mutate(Lineage2 = map(Lineage,shorten_NA))
# end_time2 <- Sys.time()
# end_time2 - start_time2
#
# taxID2lins_final <- taxID2lins_noNA %>% unite(Lineage.final, Lineage2, X3)
#
# taxID2lins_final <- taxID2lins_final %>% rename(taxid = X1)
# ### TESTING ###
# samp <- taxID2lins_collapsed[1:4,]
#
# samp_group <- rep(1:cl, length.out = nrow(samp))
#
# samp <- bind_cols(data.frame(samp_group), samp)
#
# cluster=new_cluster(cl); cluster_copy(cluster,c("map","shorten_NA", "str_locate","is.na","substr"))
#
# party = samp %>% group_by(samp_group) %>% partition(cluster)
#
# samp_noNA <- party %>% mutate(Lineage2 =shorten_NA(Lineage)) %>% collect()
#
# # samp_nNA <- samp %>% mutate(Lineage2 = map(Lineage, shorten_NA) )
#
# ####

