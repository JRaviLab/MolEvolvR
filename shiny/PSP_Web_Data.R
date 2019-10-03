library(readr)
source("R/plotting.R")
source("R/cleanup.R")
source("R/summarize.R")
#source("R/reverse_operons.R")
#source("shiny/shinyfunctions.R")

####Data import####

all <- read_tsv("data/rawdata_tsv/all_half_cleaned.txt")
# colnames(all)[colnames(all) == "GenContext"] = "GenContext.norev"
# colnames(all)[colnames(all) == "GenContext_1"] = "GenContext"
# write_tsv(all, "data/rawdata_tsv/all_half_cleaned.txt")

lineages_map <- read_delim("data/acc_files/organisms2lineages_map_bae_20170828.txt",
                           delim="\t", col_names=T, comment="#", trim_ws=T)

#domains.replace <- read_delim("data/acc_files/domains.replace.txt",
#                              delim = "\t", col_names = TRUE)

#domains.remove <- read_delim("data/acc_files/domains.ignore.txt",
#                             delim = "\t", col_names = TRUE)
####Cleanup####
#names(all)[names(all)=='GenContext'] <- "GenContext"
#all <- all %>%
  #cls_cleanup() %>%
#  cleanup_species() %>%
#  reverse_operons() %>%
  #Fix this to work on a few columns, not all
#  replace_doms_all(domains.replace,domains.remove)

all<-map(all,function(x) x %>%  		str_replace_all("\\+", " ") %>%
                     str_replace_all("(?i)\\b([a-z0-9_-]+)\\b(?:\\s+\\1\\b)+", "\\1(s)") %>%
                     str_replace_all(" ", "+")) %>% as.data.frame()

# #Rename DomArch to DomArch because that's what replace_doms requires as a column
# colnames(all)[colnames(all)=="DomArch"] <- "DomArch"
# all <- replace_doms(all,domains.replace,domains.remove)
# #Give original name back
# #all<- as.data.frame(all)
# colnames(all)[colnames(all)=="DomArch"] <- "DomArch"

####WordCounts####
create_DA.doms<- function(prot){
  DA <- prot$DomArch
  #Replace '+' with " "
  DA <- as.character(map(1:length(DA),function(x){gsub("\\+"," ", DA[x])}))
  prot <- mutate(prot,DA.doms=(DA))
  return(prot)
}
create_GC.DA <- function(prot){
  GC <- prot$GenContext
  GC <- as.character(map(1:length(GC),function(x){gsub("->|<-|\\|"," ",GC[x])}))
  GC <- str_replace_all(GC,"  +"," ")
  prot <- mutate(prot,GC.DA=GC)
  return(prot)
}



all <- create_DA.doms(all) %>%
  create_GC.DA()



DUF1700 <- all %>% filter(Query=="DUF1700-alpha-helical")
DUF1707 <- all %>% filter(Query=="DUF1707-SHOCT")
pspa <- all%>% filter(Query=="pspa")
pspb <- all%>% filter(Query=="pspb")
pspc <- all%>% filter(Query=="pspc")

pspm_table <- read_tsv("data/201905_data/pspm.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
pspm<- pspm_table%>%
  select(AccNum, Species, Lineage,
         DomArch=SIG.TM.LADB,# GenContext=GenContext,
         Length, GI, GenName, Annotation)

pspn <- all%>% filter(Query=="pspn")
liai_liaf = all%>% filter(Query=="LiaI-LiaF-TM")
toast_rack = all%>%filter(Query=="Toast-rack")

liag_table <- read_tsv("data/201905_data/liag.txt")
liag_data <- liag_table %>%reverse_operons() %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext,
         Length, GI, GenName, Annotation)

####Further steps of the process, might have to be done post filtering
# DA.doms.wc <- prot$DA.doms %>%
#   words2wc()
# GC.DA.wc <- prot$GC.DA %>%
#   words2wc()


##Also might need to be done post filter
####Lineage Summaries####
# prot.DA.summ.byLin <- summ.DA.byLin(prot)
# prot.DA.summ <- summ.DA(prot.DA.summ.byLin)
#
# ## Main Genomic Contexts -- Summarized by DA & Lineage
# prot.GC.summ.byDALin <- summ.GC.byDALin(prot)
# prot.GC.summ.byLin <- summ.GC.byLin(prot)
# prot.GC.summ <- summ.GC(prot.GC.summ.byDALin)


####Lineage plots#####
#Requires total counts
  #Total counts requires summby lin
create.DA.cummulative <- function(prot){
  colnames(prot)[colnames(prot)=="DomArch"] <- "DomArch.norep"
  colnames(prot)[colnames(prot)=="Lineage"] <- "Lineage.final"
  prot_DA.cummulative <-  summ.DA.byLin(prot) %>% total_counts( cutoff =0, type = "DA")
  colnames(prot)[colnames(prot)=="Lineage.final"] <- "Lineage"
  colnames(prot)[colnames(prot)=="DomArch.norep"] <- "DomArch"
  return(prot_DA.cummulative)
}
create.GC.cummulative <- function(prot){
  colnames(prot)[colnames(prot)=="GenContext"] <- "GenContext.norep"
  colnames(prot)[colnames(prot)=="DomArch"] <- "DomArch.norep"
  colnames(prot)[colnames(prot)=="Lineage"] <- "Lineage.final"
  prot_GC.cummulative <-  summ.GC.byDALin(prot) %>% total_counts( cutoff =0, type = "GC")
  colnames(prot)[colnames(prot)=="Lineage.final"] <- "Lineage"
  colnames(prot)[colnames(prot)=="DomArch.norep"] <- "DomArch"
  colnames(prot)[colnames(prot)=="GenContext.norep"] <- "GenContext"
  return(prot_GC.cummulative)
}
create.DA.lin <- function(prot){
  colnames(prot)[colnames(prot)=="DomArch"] <- "DomArch.norep"
  colnames(prot)[colnames(prot)=="Lineage"] <- "Lineage.final"
  prot_DA_lin <-  summ.DA.byLin(prot)
  colnames(prot)[colnames(prot)=="Lineage.final"] <- "Lineage"
  colnames(prot)[colnames(prot)=="DomArch.norep"] <- "DomArch"
  return(prot_DA_lin)
}
create.GC.lin <- function(prot){
  colnames(prot)[colnames(prot)=="DomArch"] <- "DomArch.norep"
  colnames(prot)[colnames(prot)=="Lineage"] <- "Lineage.final"
  colnames(prot)[colnames(prot)=="GenContext"] <- "GenContext.norep"
  prot_GC_lin <-  summ.GC.byLin(prot)
  colnames(prot)[colnames(prot)=="Lineage.final"] <- "Lineage"
  colnames(prot)[colnames(prot)=="GenContext.norep"] <- "GenContext"
  colnames(prot)[colnames(prot)=="DomArch.norep"] <- "DomArch"
  return(prot_GC_lin)
}
pspa_DA.cummulative <- create.DA.cummulative(pspa)
pspb_DA.cummulative <- create.DA.cummulative(pspb)
pspc_DA.cummulative <- create.DA.cummulative(pspc)
pspn_DA.cummulative <- create.DA.cummulative(pspn)
liai_liaf_DA.cummulative <- create.DA.cummulative(liai_liaf)

toast_rack_DA.cummulative <- create.DA.cummulative(toast_rack)

pspa_GC.cummulative <- create.GC.cummulative(pspa)
pspb_GC.cummulative <- create.GC.cummulative(pspb)
pspc_GC.cummulative <- create.GC.cummulative(pspc)
pspn_GC.cummulative <- create.GC.cummulative(pspn)
liai_liaf_GC.cummulative <- create.GC.cummulative(liai_liaf)

toast_rack_GC.cummulative <- create.GC.cummulative(toast_rack)

pspa_DA_lin <- create.DA.lin(pspa)
pspb_DA_lin <- create.DA.lin(pspb)
pspc_DA_lin <- create.DA.lin(pspc)
pspn_DA_lin <- create.DA.lin(pspn)
liai_liaf_DA_lin <- create.DA.lin(liai_liaf)
#liag_DA_lin <- create.DA.lin(liag)
toast_rack_DA_lin <- create.DA.lin(toast_rack)


pspa_GC_lin <- create.GC.lin(pspa)
pspb_GC_lin <- create.GC.lin(pspb)
pspc_GC_lin <- create.GC.lin(pspc)
pspn_GC_lin <- create.GC.lin(pspn)
liai_liaf_GC_lin <- create.GC.lin(liai_liaf)
toast_rack_GC_lin <- create.GC.lin(toast_rack)
# ####Upset Plot####
# #Requires DA.doms.wc or GC.DA.wc separate file
# #upset.plot(prot, 10, "da2doms")
# upset.plot(DUF1700,10,"da2doms")
# #upset.plot(prot, 75, "gc2da")
#
# colNam <-  c("AccNum", "Nghbrhd", "DomArch", "LADB","Pfam", "Length","GenName", "Lineage", "Organism", "Annotation", "GI")
#
# #DUF1700.allfa.op_ins_cls.20170817.txt
# DUF1700_table <- read_tsv(file = "data/dufs/DUF1700.allfa.op_ins_cls.20170817.txt", col_names = FALSE, skip = 1)
# colnames(DUF1700_table) <- colNam
# DUF1700_data <-  DUF1700_table%>%
#   select(AccNum, Species = Organism,
#          Lineage,
#          DomArch,# GenContext=GenContext,
#          Length, GI, GenName, Annotation)
# #DUF1707-SHOCT.1.op_ins_cls.20170817.txt
# DUF1707_table <- read_tsv(file ="data/dufs/DUF1707-SHOCT.1.op_ins_cls.20170817.txt", col_names = FALSE, skip = 1)
# colnames(DUF1707_table) <- colNam
#
# #pspa_snf7.txt
# pspa_table <- read_tsv("data/pspa_snf7.txt")
# pspa_table <- reverse_operons(pspa_table)
# pspa_data <- pspa_table %>%
#   select(AccNum, Species, Lineage=Lineage.final,
#          DomArch=DomArch, GenContext,
#          Length, GI, GenName, Annotation)
 #pspm.renamed.arch.cls.20170619.txt
# #pspb.txt
# pspb_table <- read_tsv("data/pspb.txt")
# pspb_data <- pspb_table %>% reverse_operons() %>%
#               select(AccNum, Species, Lineage=Lineage.final,
#                 DomArch=DomArch, GenContext,
#                 Length, GI, GenName, Annotation)
# #pspc.txt
# pspc_table <- read_tsv("data/pspc.txt")
# pspc_data <- pspc_table %>% reverse_operons() %>%
#   select(AccNum, Species, Lineage=Lineage.final,
#          DomArch=DomArch, GenContext,
#          Length, GI, GenName, Annotation)
# #pspn.renamed.arch.cls.20170619.txt
# pspn_table <- read_tsv("data/pspn.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
# pspn_data <- pspn_table%>%
#   select(AccNum, Species, Lineage,
#          DomArch=SIG.TM.LADB,# GenContext=GenContext,
#          Length, GI, GenName, Annotation)
# #liai.txt
# liai_table <- read_tsv("data/liai.txt")
# liai_data <- liai_table %>%reverse_operons() %>%
#   select(AccNum, Species, Lineage=Lineage.final,
#          DomArch=DomArch, GenContext,
#          Length, GI, GenName, Annotation)
# #liaf.txt
# liaf_table <- read_tsv("data/liaf.txt")
# liaf_data <- liaf_table %>%reverse_operons() %>%
#   select(AccNum, Species, Lineage=Lineage.final,
#          DomArch=DomArch, GenContext,
#          Length, GI, GenName, Annotation)
# #liag.txt

# #pspa-da_lin_counts.txt
# pspa_DA_Lin <- read_tsv("data/pspa-da_lin_counts.txt")
# #pspb-da_lin_counts
# pspb_DA_Lin <- read_tsv("data/pspb-da_lin_counts.txt")
# #pspc-da_lin_counts
# pspc_DA_Lin <- read_tsv("data/pspc-da_lin_counts.txt")
#
# liaf_DA_lin <- read_tsv("data/liaf-da_lin_counts.txt")
# liag_DA_lin <- read_tsv("data/liag-da_lin_counts.txt")
# liai_DA_lin <- read_tsv("data/liai-da_lin_counts.txt")
#
# #pspa-gc_lin_counts.txt
# pspa_GC_Lin <- read_tsv("data/pspa-gc_lin_counts.txt")
# #pspb-gc_lin_counts.txt
# pspb_GC_Lin <- read_tsv("data/pspb-gc_lin_counts.txt")
# #pspc-gc_lin_counts.txt
# pspc_GC_Lin <- read_tsv("data/pspc-gc_lin_counts.txt")
#
# liaf_GC_lin <- read_tsv("data/liaf-gc_lin_counts.txt")
# liag_GC_lin <- read_tsv("data/liag-gc_lin_counts.txt")
# liai_GC_lin <- read_tsv("data/liai-gc_lin_counts.txt")
#
# pspa_cum <- total_counts(pspa_GC_Lin,0)
# #%>%              select(GenContext, Lineage = Lineage.final, count, totalcount)
# pspb_cum <- total_counts(pspb_GC_Lin,0)
# pspc_cum <- total_counts(pspc_GC_Lin,0)
# liaf_cum <- total_counts(liaf_GC_lin,0)
# liag_cum <- total_counts(liag_GC_lin,0)
# liai_cum <- total_counts(liai_GC_lin,0)
#
# #Domain Architecture lineage with total count column
# pspa_totalC <- total_counts(pspa_DA_Lin,0, "DA")
# pspb_totalC <- total_counts(pspb_DA_Lin,0, "DA")
# pspc_totalC <- total_counts(pspc_DA_Lin,0, "DA")
# liaf_totalC <- total_counts(liaf_DA_lin,0, "DA")
# liag_totalC <- total_counts(liag_DA_lin,0, "DA")
# liai_totalC <- total_counts(liai_DA_lin,0, "DA")
# #pspa.GC.summ.byLin.v2.txt
# pspa_GC_summ_Lin <- read_tsv("data/pspa.GC.summ.byLin.v2.txt")
# #pspa.sub.v2.txt
# pspa_sub <- read_tsv("data/pspa.sub.v2.txt")
#
#
#
# pspa.DA.doms.wc <- read_tsv("data/pspa.queryDA.domains.wordcounts.v2.txt")
# pspb.DA.doms.wc <- read_tsv("data/pspb.queryDA.domains.wordcounts.txt")
# pspc.DA.doms.wc <- read_tsv("data/pspc.queryDA.domains.wordcounts.v1.txt")
# liaf.DA.doms.wc <- read_tsv("data/liaf.queryDA.domains.wordcounts.txt")
# liag.DA.doms.wc <- read_tsv("data/liag.queryDA.domains.wordcounts.txt")
# liai.DA.doms.wc <- read_tsv("data/liai.queryDA.domains.wordcounts.txt")
#
# pspa.GC.doms.wc <- read_tsv("data/pspa.GC.DA.wordcounts.v2.txt")
# pspb.GC.doms.wc <- read_tsv("data/pspb.GC.DA.wordcounts.txt")
# pspc.GC.doms.wc <- read_tsv("data/pspc.GC.DA.wordcounts.v1.txt")
# liaf.GC.doms.wc <- read_tsv("data/liaf.GC.DA.wordcounts.txt")
# liag.GC.doms.wc <- read_tsv("data/liag.GC.DA.wordcounts.txt")
# liai.GC.doms.wc <- read_tsv("data/liai.GC.DA.wordcounts.txt")
