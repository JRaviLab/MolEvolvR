## Functions to convert cluster files .op_ins_cls --> tsv
## Created: Aug 08, 2019
## Modified: Dec 11, 2019 (chunk moved from the_process.Rmd)
## Janani Ravi (@jananiravi) & Samuel Chen (@samuelzornchen)

## Run this script only if you do not have:
## 'all.op_ins_cls.clus2table' OR 'all.txt'
## Also, useful if you have individual files for each protein/domain.

##!! Note: this may not work at all since the column names in
## op_ins_cls & clean_clust_file.R may not match. Check files to fix issues!

#################
## Pkgs needed ##
#################
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")

###############################
## Convert op_ins_cls to tsv ##
###############################
## Initialize dataframe for the combined dataset; Will serve for even 1.
all <- data.frame(matrix(ncol=11, nrow=0))
colnames(all) <- colnames.op_ins_cls

for(x in list.files("data/rawdata_opinscls")){
  print(x)
  inpath <- paste0("data/rawdata_opinscls/", x)
  prot_name <- str_remove(x, ".op_ins_cls")
  # clean_clust_file(path, writepath=NULL, query)
  prot_opinscls <- clean_clust_file(path=inpath, query=prot_name)

  all <- bind_rows(all, prot_opinscls)
  # comment next lines to avoid overwriting each time!
  # Last written: Aug 20, 2019
  # write_tsv(prot_opinscls,
  #           path=paste0("data/rawdata_tsv/", prot_name, ".txt"),
  #           col_names=T)
}

# Appending genome-specific signature columns.
# Reading TaxIDs
all_taxid <- read_tsv("data/acc_files/all.taxid", col_names=F)
colnames(all_taxid) <- c("AccNum", "TaxID", "Species.q")
# Reading GCA numbers
all_gca <- read_tsv("data/acc_files/all.gca", col_names=F)
colnames(all_gca) <- c("AccNum", "GCA_ID")
# Adding TaxIDs and GCA_IDs to combined "all" dataframe
all <- all %>%
  left_join(all_taxid, by="AccNum") %>%
  left_join(all_gca, by="AccNum")
# Rearrange columns for easy readability
all <- all %>%
  select(AccNum, Query, ClustID, ClustName,
         DomArch.orig, GenContext.orig,
         Lineage, Species.q, TaxID, GCA_ID,
         arch.PFAM, arch.TMSIG,
         Length, GenName,
         Species.orig, Annotation, GI)

# Write the new combined file; Comment next line(s) to avoid overwriting each time!
# Last written: Aug 20, 2019
# write_tsv(all, path="data/rawdata_tsv/all_raw.txt", col_names=T)
