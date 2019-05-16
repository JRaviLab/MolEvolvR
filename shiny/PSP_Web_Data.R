library(readr)
source("shinyfunctions.R")

colNam <-  c("AccNum", "Nghbrhd", "DomArch", "LADB","Pfam", "Len","GenNam", "Lineage", "Organism", "Annotation", "GI")

#DUF1700.allfa.op_ins_cls.20170817.txt
DUF1700_table <- read_tsv(file = "data/dufs/DUF1700.allfa.op_ins_cls.20170817.txt", col_names = FALSE, skip = 1)
colnames(DUF1700_table) <- colNam
#DUF1707-SHOCT.1.op_ins_cls.20170817.txt
DUF1707_table <- read_tsv(file ="data/dufs/DUF1707-SHOCT.1.op_ins_cls.20170817.txt", col_names = FALSE, skip = 1)
colnames(DUF1707_table) <- colNam

#pspa_snf7.txt
pspa_table <- read_tsv("data/pspa_snf7.txt")
pspa_data <- pspa_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#pspm.renamed.arch.cls.20170619.txt
pspm_table <- read_tsv("data/pspm.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
pspm_data <- pspm_table%>%
  select(AccNum, Species, Lineage,
         SIG.TM.LADB,# GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#pspb.txt
pspb_table <- read_tsv("data/pspb.txt")
pspb_data <- pspb_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#pspc.txt
pspc_table <- read_tsv("data/pspc.txt")
pspc_data <- pspc_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#pspn.renamed.arch.cls.20170619.txt
pspn_table <- read_tsv("data/pspn.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
pspn_data <- pspn_table%>%
  select(AccNum, Species, Lineage,
         SIG.TM.LADB,# GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#liai.txt
liai_table <- read_tsv("data/liai.txt")
liai_data <- liai_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#liaf.txt
liaf_table <- read_tsv("data/liaf.txt")
liaf_data <- liaf_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)
#liag.txt
liag_table <- read_tsv("data/liag.txt")
liag_data <- liag_table %>%
  select(AccNum, Species, Lineage=Lineage.final,
         DomArch=DomArch.norep, GenContext=GenContext.norep,
         Length, GI, GenName, Annotation)

#pspa-da_lin_counts.txt
pspa_DA_Lin <- read_tsv("data/pspa-da_lin_counts.txt")
#pspb-da_lin_counts
pspb_DA_Lin <- read_tsv("data/pspb-da_lin_counts.txt")
#pspc-da_lin_counts
pspc_DA_Lin <- read_tsv("data/pspc-da_lin_counts.txt")

liaf_DA_lin <- read_tsv("data/liaf-da_lin_counts.txt")
liag_DA_lin <- read_tsv("data/liag-da_lin_counts.txt")
liai_DA_lin <- read_tsv("data/liai-da_lin_counts.txt")

#pspa-gc_lin_counts.txt
pspa_GC_Lin <- read_tsv("data/pspa-gc_lin_counts.txt")
#pspb-gc_lin_counts.txt
pspb_GC_Lin <- read_tsv("data/pspb-gc_lin_counts.txt")
#pspc-gc_lin_counts.txt
pspc_GC_Lin <- read_tsv("data/pspc-gc_lin_counts.txt")

liaf_GC_lin <- read_tsv("data/liaf-gc_lin_counts.txt")
liag_GC_lin <- read_tsv("data/liag-gc_lin_counts.txt")
liai_GC_lin <- read_tsv("data/liai-gc_lin_counts.txt")

pspa_cum <- cummulative.count(pspa_GC_Lin,0)
pspb_cum <- cummulative.count(pspb_GC_Lin,0)
pspc_cum <- cummulative.count(pspc_GC_Lin,0)
liaf_cum <- cummulative.count(liaf_GC_lin,0)
liag_cum <- cummulative.count(liag_GC_lin,0)
liai_cum <- cummulative.count(liai_GC_lin,0)

#Domain Architecture lineage with total count column
pspa_totalC <- cummulative.count(pspa_DA_Lin,0, "DA")
pspb_totalC <- cummulative.count(pspb_DA_Lin,0, "DA")
pspc_totalC <- cummulative.count(pspc_DA_Lin,0, "DA")
liaf_totalC <- cummulative.count(liaf_DA_lin,0, "DA")
liag_totalC <- cummulative.count(liag_DA_lin,0, "DA")
liai_totalC <- cummulative.count(liai_DA_lin,0, "DA")
#pspa.GC.summ.byLin.v2.txt
pspa_GC_summ_Lin <- read_tsv("data/pspa.GC.summ.byLin.v2.txt")
#pspa.sub.v2.txt
pspa_sub <- read_tsv("data/pspa.sub.v2.txt")



pspa.DA.doms.wc <- read_tsv("data/pspa.queryDA.domains.wordcounts.v2.txt")
pspb.DA.doms.wc <- read_tsv("data/pspb.queryDA.domains.wordcounts.txt")
pspc.DA.doms.wc <- read_tsv("data/pspc.queryDA.domains.wordcounts.v1.txt")
liaf.DA.doms.wc <- read_tsv("data/liaf.queryDA.domains.wordcounts.txt")
liag.DA.doms.wc <- read_tsv("data/liag.queryDA.domains.wordcounts.txt")
liai.DA.doms.wc <- read_tsv("data/liai.queryDA.domains.wordcounts.txt")

pspa.GC.doms.wc <- read_tsv("data/pspa.GC.DA.wordcounts.v2.txt")
pspb.GC.doms.wc <- read_tsv("data/pspb.GC.DA.wordcounts.txt")
pspc.GC.doms.wc <- read_tsv("data/pspc.GC.DA.wordcounts.v1.txt")
liaf.GC.doms.wc <- read_tsv("data/liaf.GC.DA.wordcounts.txt")
liag.GC.doms.wc <- read_tsv("data/liag.GC.DA.wordcounts.txt")
liai.GC.doms.wc <- read_tsv("data/liai.GC.DA.wordcounts.txt")
