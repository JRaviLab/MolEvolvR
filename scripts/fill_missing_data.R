#####Add Missing data Samuel Chen #####

setwd("C:/Users/samue/Google_Drive/GitHub/the-approach")
source("R/cleanup.R")
# source("R/cleanup_domarch.R")
source("R/clean_clust_file.R")
# Initialize dataframe for the combined dataset; Will serve for even 1.
all <- data.frame(matrix(ncol=11, nrow=0))
colnames(all) <- colnames.op_ins_cls

for(x in list.files("data/rawdata_opinscls")){
  print(x)
  inpath <- paste0("data/rawdata_opinscls/", x)
  prot_name <- str_remove(x, ".op_ins_cls")
  prot_opinscls <- clean_clust_file(path=inpath, query=prot_name)

  all <- bind_rows(all, prot_opinscls)
}
all <- filter(all, AccNum!="-")

# Reading TaxIDs
all_taxid <- read_tsv("data/acc_files/all.taxid", col_names=F)
colnames(all_taxid) <- c("AccNum", "TaxID", "Species.q")

#Remove the NA duplicates of distinct_tax
reduced_all_taxid <- filter(all_taxid, !(AccNum=="ABX48815.1"&is.na(TaxID))
                            & !(AccNum=="AJK88848.1" &is.na(TaxID))
                            & !(AccNum=="AKK03387.1" &is.na(TaxID)) ) %>% distinct()

# Reading GCA numbers
all_gca <- read_tsv("data/acc_files/all.gca", col_names=F)
colnames(all_gca) <- c("AccNum", "GCA_ID")

#This is the ideal, no added rows
all_reduce_merge <- all %>%
  left_join(reduced_all_taxid, by="AccNum")

all_reduce_merge <- all_reduce_merge%>%
  left_join(distinct(all_gca),by="AccNum")


####Use df below to fill the rows that are missing an observations
missing_gis <-read_tsv("data/acc_files/prob_files/missing_data_ba.txt", col_names = F)
colnames(missing_gis) = c("AccNum", "Length", "GeneName", "Lineage", "Species", "GCA_ID", "Annotation","GI")

#How do I match the rows? AccNum column?
##Compile list of AccNum from missing_gis, then
##Use AccNum to filter, if there are more than one AccNum, go to
all_mis <- missing_gis %>% group_by(AccNum) %>% count() %>% arrange(-n) # AccNum unique in missing
all_acc <- all_reduce_merge %>%group_by(AccNum) %>%count() %>% arrange(-n) %>% filter(n!=1)


missing_spp_taxids <-read_tsv("data/acc_files/prob_files/missing_spp_taxids.txt", col_names = F)
colnames(missing_spp_taxids) = c("AccNum", "TaxID", "Species")
missing_spp_taxids %>% group_by(AccNum,Species) %>% count() %>% arrange(-n)

missing_gis <- left_join(missing_gis,missing_spp_taxids, c("AccNum","Species") ) %>% distinct()#11 extra rows

length(setdiff(all_mis$AccNum,all_acc$AccNum)) #only 454, I want it to be 457
#intersect(all_mis$AccNum,all_acc$AccNum)
#ABX48815.1 is the same in both instances of all, other than Query and clustid
#### Why is gen_context missing from one of them?
##Many of the missing values are not filled in by table/ have no effect on respective rows
all_test <- all_reduce_merge
for( x in missing_gis$AccNum){
  acc_match <- filter(missing_gis, AccNum==x) #Gets the observations of row
  acc_index <- grep(x,all_reduce_merge$AccNum) #Gets the row numbers
  all_test[acc_index,"Length"] <- acc_match$Length
  all_test[acc_index,"GenName"] <- acc_match$GeneName
  all_test[acc_index,"Lineage"] <- acc_match$Lineage
  all_test[acc_index,"Species.q"] <- acc_match$Species
  all_test[acc_index,"GCA_ID"] <- acc_match$GCA_ID
  all_test[acc_index,"Annotation"] <- acc_match$Annotation
  all_test[acc_index,"GI"] <- acc_match$GI
  all_test[acc_index,"TaxID"] <- acc_match$TaxID
}


diff <- setdiff(all_test,all_reduce_merge)

oddballs <- missing_gis %>% group_by(Lineage) %>% count()

##Find the duplicates in same Query
all_reduce_merge %>% group_by(AccNum, Query) %>% count() %>% arrange(-n) #all ones so we good

missing_gis_gca_spp_taxids <-read_tsv("data/acc_files/prob_files/missing_gis_gca_spp_taxids_foraccnum.txt")

#write_tsv(all_reduce_merge, "data/rawdata_tsv/all_merged_gca_taxid.txt",col_names = T)
