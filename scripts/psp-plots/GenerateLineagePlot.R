library(tidyverse)
source("R/plotting.R")


all <- read_tsv("data/rawdata_tsv/all_clean.txt")

domains_of_interest <- c("ZNR", "SIGMA-HTH", "GNTR-HTH", "Beta-Propeller", "DUF2089-HTH", "SHOCT-like",
                         "DUF1707−SHOCT−bihelical", # This is currently just DUF1707-SHOCT: rename
                         "MacB_PCD", "FTSW_RODA_SPOVE", "PADR−HTH", "DUF1700−alpha−helical",
                         "REC", "HISKIN", "LiaI−LiaF−TM", "Toast-rack", "PspC", "PspB", "PspF",
                         "Tfu_10090", "Psp-AA", "Spermine_synth", "Flot", "Band_7", "Cest_Tir",
                         "PspA", "Snf7", "Classical-AAA", "TraJ-RHH", "clgR-HTH", "Thioredoxin",
                         "DUF3046","PspN_N", "PspM")


all$DomArch <- purrr::map(all$DomArch, function(x)str_replace_all(x, pattern = "DUF1707-SHOCT", "DUF1707-SHOCT-bihelical" ))

all_grouped <- data.frame("Query" = character(0), "Lineage" = character(0), "count"= integer())
for(dom in domains_of_interest)
{
  domSub <- all %>% filter(grepl(dom, DomArch, ignore.case = T))

  domSub <- domSub %>% group_by(Lineage) %>% summarize("count" = n())

  domSub$Query = dom

  all_grouped <- dplyr::union(all_grouped, domSub)
}

GetKingdom <- function(lin)
{
  gt_loc <- str_locate(lin, ">")[,"start"]

  if(is.na(gt_loc))
  {
    # No '>' in lineage
    return(lin)
  }
  else
  {
    kingdom <- substring(lin, first = 0, last = (gt_loc-1))
    return(kingdom)
  }
}

GetPhylum <- function(lin)
{
  gt_loc <- str_locate_all(lin, ">")[[1]]

  if(length(gt_loc) == 0)
  {
    # No '>' in lineage
    return(lin)
  }
  else
  {
    gt_loc <- gt_loc[,"start"][length(gt_loc[,"start"])]
    kingdom <- substring(lin, first = (gt_loc+1) )
    return(kingdom)
  }
}

all_grouped <- all_grouped %>% mutate(Kingdom = unlist(purrr::map(Lineage,GetKingdom)) ) %>%
  mutate(Phylum = unlist(purrr::map(Lineage, GetPhylum)))

all_grouped$Lineage <- factor(all_grouped$Lineage,
                              levels = sort(names(sort(table(all_grouped$Lineage),decreasing=TRUE))))


lin_counts <- all_grouped %>% group_by(Kingdom, Phylum) %>% summarize("count" =n())

grep("bacteria", lin_counts$Kingdom) %>% length()

bacteria_colors <- rep("#d94f25", length(grep("bacteria", lin_counts$Kingdom)))

archaea_colors <- rep("#26662d", length(grep("archaea", lin_counts$Kingdom)))

eukaryota_colors <- rep("#123d99", length(grep("eukaryota", lin_counts$Kingdom)))

virus_colors <- rep("#538073", length(grep("virus", lin_counts$Kingdom)))
colors <- append(archaea_colors, bacteria_colors) %>%
  append(eukaryota_colors) %>% append(virus_colors)


# all_grouped <- all_grouped %>% mutate(Lineage = unlist(purrr::map(Lineage, GetPhylum)))

ggplot(data=all_grouped,
       aes_string(x="Lineage", y="Query")) +
  geom_tile(data=subset(all_grouped,
                        !is.na(count)),
            aes(fill=count),
            colour="darkred", size=0.3) + #, width=0.7, height=0.7),
  scale_fill_gradient(low="white", high="darkred") +
  # scale_x_discrete(position="top") +
  theme_minimal()+
  theme(axis.title = element_blank() ,axis.text.x=element_text(angle=90,hjust=1,vjust=0.2, color =
                                   colors
                                   ))




