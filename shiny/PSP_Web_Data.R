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
#pspm.renamed.arch.cls.20170619.txt
pspm_table <- read_tsv("data/pspm/pspm.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
#pspb.txt
pspb_table <- read_tsv("data/pspb.txt")
#pspc.txt
pspc_table <- read_tsv("data/pspc.txt")
#pspn.renamed.arch.cls.20170619.txt
pspn_table <- read_tsv("data/pspn-duf3046/pspn.renamed.arch.cls.20170619.txt") %>% filter(!grepl(pattern =  "#", AccNum))
#liai.txt
liai_table <- read_tsv("data/liai.txt")
#liaf.txt
liaf_table <- read_tsv("data/liaf.txt")
#liag.txt
liag_table <- read_tsv("data/liag.txt")

#pspa-da_lin_counts.txt
pspa_DA_Lin <- read_tsv("data/pspa-da_lin_counts.txt")
#pspb-da_lin_counts
pspb_DA_Lin <- read_tsv("data/pspb-da_lin_counts.txt")
#pspc-da_lin_counts
pspc_DA_Lin <- read_tsv("data/pspc-da_lin_counts.txt")

#pspa-gc_lin_counts.txt
pspa_GC_Lin <- read_tsv("data/pspa-gc_lin_counts.txt")
#pspb-gc_lin_counts.txt
pspb_GC_Lin <- read_tsv("data/pspb-gc_lin_counts.txt")
#pspc-gc_lin_counts.txt
pspc_GC_Lin <- read_tsv("data/pspc-gc_lin_counts.txt")

pspa_cum <- cummulative.count(pspa_GC_Lin,0)
pspb_cum <- cummulative.count(pspb_GC_Lin,0)
pspc_cum <- cummulative.count(pspc_GC_Lin,0)

#Domain Architecture lineage with total count column
pspa_totalC <- cummulative.count(pspa_DA_Lin,0, "DA")
pspb_totalC <- cummulative.count(pspb_DA_Lin,0, "DA")
pspc_totalC <- cummulative.count(pspc_DA_Lin,0, "DA")

#pspa.GC.summ.byLin.v2.txt
pspa_GC_summ_Lin <- read_tsv("data/pspa-snf7/pspa.GC.summ.byLin.v2.txt")
#pspa.sub.v2.txt
pspa_sub <- read_tsv("data/pspa-snf7/pspa.sub.v2.txt")



pspa.DA.doms.wc <- read_tsv("data/pspa-snf7/pspa.queryDA.domains.wordcounts.v2.txt")
#pspb.doms.wc <- read_tsv("")
pspc.DA.doms.wc <- read_tsv("data/pspbc/pspc.queryDA.domains.wordcounts.v1.txt")

pspa.GC.doms.wc <- read_tsv("data/pspa-snf7/pspa.GC.DA.wordcounts.v2.txt")
#pspb.GC.doms.wc <-
pspc.GC.doms.wc <- read_tsv("data/pspbc/pspc.GC.DA.wordcounts.v1.txt")
