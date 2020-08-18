library(tidyverse)
library(rlang)
source("R/network.R")
source("R/plotting.R")

all <- read_tsv("data/rawdata_tsv/all_clean.txt")

queries = c("PspA", "Snf7", "PspB","PspC", "PspM", "PspN","DUF3046", "LiaI-LiaF-TM", "Toast-rack",
            "Tfu_1009","DUF1700", "DUF1707")


#### heatmaps ####

png(file = "data/figures/poster_figs/lineage_query.png", width = 3500, height = 2000,res = 300)

lineage.Query.plot(all, queries, "DomArch", 100)
dev.off()


png(file = "data/figures/poster_figs/lineage_DomainArchitecture.png", width = 4000, height = 2000,res = 300)
lineage.DA.plot(all, "DomArch", cutoff = 90)
dev.off()

png(file = "data/figures/poster_figs/lineage_GenContext.png", width = 4000, height = 2000,res = 300)
lineage.DA.plot(pspc, "GenContext", cutoff = 50)
dev.off()

png(file = "data/figures/poster_figs/viridis_lineage_DomainArchitecture.png", width = 4000, height = 2000,res = 300)
lineage.DA.plot(all, "DomArch", cutoff = 90, color = "viridis")
dev.off()

png(file = "data/figures/poster_figs/viridis_lineage_GenContext.png", width = 4000, height = 2000,res = 300)
lineage.DA.plot(all, "GenContext", cutoff = 30, color = "viridis")
dev.off()

pspa <- all %>% filter(grepl("PspA", DomArch,ignore.case = T))

#### Network ####
# png(file = "data/figures/poster_figs/domainNetwork.png", width = 4000, height = 2000,res = 300)
domain_network(all, "DomArch.repeats", queries, 90)
# dev.off()

#### Upset Plots ####
png(file = "data/figures/poster_figs/upsetDomainArchitecture.png", width = 4000, height = 2000,res = 300)
upset.plot(toast_rack, "DomArch", cutoff = 100)
dev.off()

png(file = "data/figures/poster_figs/upsetGenomicContext.png", width = 2000, height = 900,res = 300)
upset.plot(all, "GenContext", cutoff = 25)
dev.off()

#### Sunburst ####
# png(file = "data/figures/poster_figs/pspaSunburst.png", width = 4000, height = 2000,res = 300)
lineage_sunburst(pspa, "Lineage", levels = 2)
# dev.off()

#### Stacked Barplot ####
prot_data <- total_counts(all, column = "GenContext", cutoff = 40)

# DAorder used to order bars in descending order
GCorder <- prot_data %>% select(GenContext) %>% unique() %>% pull(GenContext) %>% rev()

prot_data <- prot_data %>% select(GenContext, Lineage, count)

prot_data$Lineage = unlist(map(prot_data$Lineage, function(x){
  gt_pos = gregexpr(pattern = ">", x)[[1]][2]
  # If there is not a second '>' in the lineage
  if(is.na(gt_pos))
  {
    x
  }
  else{
    substr(x,1, (gt_pos-1) )
  }
} ))

prot_data$GenContext <- factor(prot_data$GenContext, levels = GCorder)

png(file = "data/figures/poster_figs/stackedPlot.png", width = 5000, height = 2000,res = 300)
(ggplot(prot_data, aes(fill = Lineage, y = count, x = GenContext)) +
  geom_bar(position = 'stack', stat = "identity") +
  coord_flip() +
  xlab("Genomic Context")+
  ylab("Count") +
  theme_minimal() + theme(legend.position = c(0.7, 0.4)))
dev.off()







