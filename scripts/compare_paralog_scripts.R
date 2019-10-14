#setwd("C:/Users/samue/Google_Drive/GitHub/the-approach")

# Also, get some numbers for eukaryotes!

find_paralogs_new <- function(df){
  #Remove eukaryotes
  df <- df %>% filter(!grepl("^eukaryota",Lineage))
  paralogTable <- df %>%
    group_by(Query,GCA_ID) %>%
    count(DomArch.orig) %>%
    filter(n>1) %>%
    arrange(-n)
  colnames(paralogTable)[colnames(paralogTable)=="n"] = "Count"
  ###Merge with columns: AccNum,TaxID, and GCA/ Species?
  #paralogTable <- df %>% select(AccNum, TaxID, GCA_ID) %>%
  paralogTable <- df %>% select(Species.orig , TaxID,  GCA_ID, Query) %>%
    left_join(paralogTable, by= c("GCA_ID","Query")) %>%
    filter(!is.na(Count)) %>% distinct()
  return(paralogTable)
}

find_paralogs_old <- function(df){
  #Remove eukaryotes
  df <- df %>% filter(!grepl("^eukaryota",Lineage))
  paralogTable <- df
  print(colnames(paralogTable))
  paralogTable <-  group_by(.data = paralogTable, Query,Species.orig ) %>%
    count()%>%
    filter(n>1) %>% arrange(-n)
  colnames(paralogTable)[colnames(paralogTable)=="n"] = "Count"
  ###Merge with columns: AccNum,TaxID, and GCA/ Species?
  return(paralogTable)
}

source("R/summarize.R")

all_with_tax_gca <- read_tsv("data/rawdata_tsv/all_merged_gca_taxid.txt")
df <- all_with_tax_gca
para1 <- find_paralogs_old(all_with_tax_gca)
para2 <- find_paralogs_new(all_with_tax_gca)


para_with_euk <- find_paralogs(all_with_tax_gca)

no_gca <- filter(all_with_tax_gca, GCA_ID == "-")
euk_nogca <- filter(no_gca, grepl("^eukaryota" ,Lineage))
non_euk <- filter(no_gca, !grepl("^eukaryota" ,Lineage))


para_no_euk <- find_paralogs(all_with_tax_gca)
