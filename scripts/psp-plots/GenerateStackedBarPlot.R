library(tidyverse)

source("R/cleanup.R")

# read in all and clean
all <- read_tsv("data/rawdata_tsv/all_semiclean.txt")

colnames(all)[which(colnames(all) == "DomArch")] = "DomArch.ind"
colnames(all)[which(colnames(all) == "DomArch.orig")] = "DomArch.ind.orig"

colnames(all)[which(colnames(all) == "ClustName")] = "DomArch"

all <- all %>%
  cleanup_clust(domains_keep = domains_keep, domains_rename = domains_rename,repeat2s = FALSE)
colnames(all)[which(colnames(all) == "ClustName")] = "DomArch.repeats"

colnames(all)[which(colnames(all) == "ClustName.orig")] = "DomArch.orig"

colnames(all)[which(colnames(all) == "GenContext")] = "Temp"

all <- cleanup_gencontext(all,domains_rename = domains_rename, repeat2s = FALSE)
colnames(all)[which(colnames(all) == "GenContext")] = "GenContext.repeats"
colnames(all)[which(colnames(all) =="Temp")] = "GenContext"


all <- all %>% select(AccNum, DomArch, DomArch.repeats, GenContext, GenContext.repeats,
                      Lineage, Species, GeneName, Length, GCA_ID) %>% distinct()




## data frame requires DomArch(Clustname) and Lineage

prot <- toast_rack# all %>% filter(!grepl("PspA|Snf7", DomArch))


# x = DomArch
# fill = lineage
# y = counts of each

# DomArch  Lineage  Counts

prot_data <- total_counts(prot, column = "DomArch", cutoff = 90)

# DAorder used to order bars in descending order
DAorder <- prot_data %>% select(DomArch) %>% unique() %>% pull(DomArch) %>% rev()

prot_data <- prot_data %>% select(DomArch, Lineage, count)

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

prot_data$DomArch <- factor(prot_data$DomArch, levels = DAorder)

perc <- ggplot(prot_data, aes(fill = Lineage, y = count, x = DomArch)) +
  geom_bar(position = 'fill', stat = "identity") +
  coord_flip() +
  xlab("Domain Architecture")+
  ylab("Percent") +
  theme_minimal()
  # scale_fill_viridis(discrete = T)
stack <- ggplot(prot_data, aes(fill = Lineage, y = count, x = DomArch)) +
  geom_bar(position = 'stack', stat = "identity") +
  coord_flip() +
  xlab("Domain Architecture")+
  ylab("Counts") +
  theme_minimal() + theme(legend.position = c(0.7, 0.4))
## Increase axes and legend font size

stack
# perc
