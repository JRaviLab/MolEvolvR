## Script to find and remove duplicate entries to get accurate counts
## Old starting repo: `the_approach`
## Tested with PSP data
## Need to generalize a bit more & convert to a function

## Created: Oct 2019
## Modified: Nov 2022 | Move copy from molevol_scripts to psp_app
## Authors: Samuel Chen, Janani Ravi

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

#No difference in the AccNum columns of the 2
#setdiff(all_taxid$AccNum,all$AccNum)


#####Summarise All, Get the duplicates and the counts
all_summarise <- all %>% group_by(AccNum)
all_count <- all_summarise %>% count() %>% arrange(-n)

all_dup <- all_count %>% filter(n>=2)
all_count_vec <- data.frame(c("Count_4"=length(filter(all_dup,n==4)$n),
                              "Count_3"=length(filter(all_dup,n==3)$n),
                              "Count_2"=length(filter(all_dup,n==2)$n)))

# Reading TaxIDs
all_taxid <- read_tsv("data/acc_files/all.taxid", col_names=F)
colnames(all_taxid) <- c("AccNum", "TaxID", "Species.q")

####Merge All and the taxid and figure
all_tax <- all %>%
  left_join(all_taxid, by="AccNum")

##Merge an all without duplicates with alltax: result in 21797 obs
#all_wo_dup <- select(all,AccNum) %>% distinct()
#all_wo_dup <- mutate(all_wo_dup,ColumnHolder = row.names(all_wo_dup)) %>%
  left_join(all_taxid,by="AccNum")

##Merge an alltax without duplicates with all: result in 21797 obs
#alltax_wo_dup <- select(all_taxid,AccNum) %>% distinct()
#alltax_wo_dup <- mutate(alltax_wo_dup,ColumnHolder = row.names(alltax_wo_dup)) %>% left_join(all,by="AccNum")

#Merg
#alltax_wo_dup <- select(all_taxid,AccNum) %>% distinct()
#alltax_wo_dup <- mutate(alltax_wo_dup,ColumnHolderTax = paste0(row.names(alltax_wo_dup),"tax"))
#all_wo_dup <- select(all,AccNum) %>% distinct()
#all_wo_dup <- mutate(all_wo_dup,ColumnHolder = paste0(row.names(all_wo_dup),"All")) %>%
  left_join(alltax_wo_dup,by="AccNum")

######Remove one of each quad, triple, and double and see how rows react ######
# #Remove a 4 count and test how tax reacts
# #: 21793 -> 24369 : loses 16 total: loses 4 original all
# all_no4 <-filter(all,AccNum!="ADQ78255.1")
# all_taxid_no4 <-filter(all_taxid,AccNum!="ADQ78255.1")
# all_tax_no4 <- all_no4 %>%
#   left_join(all_taxid_no4, by="AccNum")
##Remove a 3 count and test how tax reacts
##: 21794 -> 24376 : loses 9 total: 3*3 : loses 3 original all
# all_no3 <-filter(all,AccNum!="ABQ03233.1")
# all_taxid_no3 <-filter(all_taxid,AccNum!="ABQ03233.1")
# all_tax_no3 <- all_no3 %>%
#   left_join(all_taxid_no3, by="AccNum")

##Remove a 2 count and test how tax reacts
#: 21795 -> 24381 : loses 4 total: 2*2 : loses 2 original
# all_no2 <-filter(all,AccNum!="ABJ67228.1")
# all_taxid_no2 <-filter(all_taxid,AccNum!="ABJ67228.1")
# all_tax_no2 <- all_no2 %>%
#   left_join(all_taxid_no2, by="AccNum")
##9*4+83*3+991*2 = 2267,   Total Rows gained is 2588

######Create a taxid df with 1 obs. per AccNum
distinct_tax <- all_taxid %>% distinct() %>% group_by(AccNum) %>% count() %>% arrange(-n)

#Remove the NA duplicates of distinct_tax
reduced_all_taxid <- filter(all_taxid, !(AccNum=="ABX48815.1"&is.na(TaxID))
                        & !(AccNum=="AJK88848.1" &is.na(TaxID))
                        & !(AccNum=="AKK03387.1" &is.na(TaxID)) ) %>% distinct()
#This is the ideal, no added rows
all_reduce_merge <- all %>%
  left_join(reduced_all_taxid, by="AccNum")
#Can we use this to join?


# Reading GCA numbers
all_gca <- read_tsv("data/acc_files/all.gca", col_names=F)
colnames(all_gca) <- c("AccNum", "GCA_ID")

#No duplicates after distinct
distinct_gca <- distinct(all_gca) %>% group_by(AccNum)# %>% count() %>% arrange(-n)

all_reduce_merge <- all_reduce_merge%>%
  left_join(distinct_gca,by="AccNum")

#######Find Na's in GCA and TaxID########
#Create list with NA's in the taxid columns
#figure out which ones of those have other AccNum obs.

taxid_na <- all_taxid[which(is.na(all_taxid$TaxID)| is.na(all_taxid$Species.q) | is.na(all_taxid$AccNum)),]

No_NA <-anti_join(all_taxid, taxid_na)
matches <- paste(taxid_na$AccNum,collapse="|")

obs_and_na <- No_NA[grep(matches, No_NA$AccNum),]

###GCA NA###
##No NA's, use "-"
#gca_na <- all_gca[which(is.na(all_gca$GCA_ID)| is.na(all_gca$AccNum)),]
gca_na <- all_gca %>% filter(grepl("^-$",all_gca$GCA_ID))
