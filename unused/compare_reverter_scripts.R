## Script to compare reverter functions for genomic contexts (GenContext)
## Old starting repo: `the_approach`
## Tested with PSP data

## Created: Oct 2019
## Authors: Samuel Chen, Janani Ravi

source("R/reverse_operons.R")
source("scripts/reverse_operons.v1.R")
source("R/clean_clust_file.R")


#Compare original summarize vs reverter 2 summarize
#apply summary function to all 3 and push into github
source("R/summarize.R")
all<- read_tsv("data/rawdata_tsv/all_merged_gca_taxid.txt")
#pspa_test <- select(query, GenContext.norep)
#colnames(pspa_test) <- c("GenContext.orig")
#pspa_testrev1 <- reverse_operons(pspa_test)

all_rev1 <- select(all,GenContext.orig) %>% reverse_operons()
all_rev1_summarize <- all_rev1 %>% group_by(GenContext.norep) %>% count() %>% arrange(-n)

all_rev2 <- select(all,GenContext.orig) %>% reverse_operon(all$GenContext.orig)
all_rev2_summarize <- all_rev2 %>% group_by(GenContext.norep) %>% count() %>% arrange(-n)

#write_tsv(all_origin,"C:/Users/samue/Google_Drive/GitHub/the-approach/data/acc_files/prob_files/all_gc_origi.txt")
#write_tsv(all_rev1_summarize,"C:/Users/samue/Google_Drive/GitHub/the-approach/data/acc_files/prob_files/all_gc_rev1.txt")
#write_tsv(all_rev2_summarize,"C:/Users/samue/Google_Drive/GitHub/the-approach/data/acc_files/prob_files/all_gc_rev2.txt")

rev1 <- read_tsv("data/acc_files/prob_files/all_gc_rev1.txt")
rev2 <- read_tsv("data/acc_files/prob_files/all_gc_rev2.txt")

d1 <- setdiff(rev1$GenContext.norep,rev2$GenContext.norep)

setdiff(rev2$GenContext.norep,rev1$GenContext.norep)

anti_rev <- anti_join(rev1,rev2, "GenContext.norep")

anti_rev2 <- anti_join(rev2,rev1, "GenContext.norep")

#write_tsv(anti_rev, "C:/Users/samue/Google_Drive/GitHub/the-approach/data/acc_files/prob_files/rev1_antijoin.txt")
#write_tsv(anti_rev2, "C:/Users/samue/Google_Drive/GitHub/the-approach/data/acc_files/prob_files/rev2_antijoin.txt")


all_na_spec <- all %>% filter(is.na(Species.q))# & is.na(Species.orig))
all_na_spec <- all %>% filter(is.na(Species.orig))
