library(tidyverse)
source("R/plotting.R")

# all <- read_tsv("data/rawdata_tsv/all_clean.txt")
#
# # DUF1700-ahelical -> ahelical
# # DUF1700-like_alpha_helical_domain	 -> DUF1700-ahelical
# # DUF1707-SHOCT-bihelical ->	SHOCT-bihelical
# # Toast-rack -> Toast_rack
# # PADR-HTH	PadR-HTH
# # clgR-HTH	ClgR-HTH
# # ZNR	ZnR
# # DUF1700-like_alpha_helical_domain	DUF1700-ahelical
# # Psp-AA	PspAA
# # pspF	PspF
# # Flot	TM-Flotillin
# # Phageshock_PspD	PspD
#
# all$GenContext.repeats <- all$GenContext.repeats %>% str_replace_all("DUF1700-ahelical", "ahelical") %>%
#                       str_replace_all("DUF1700-like_alpha_helical_domain", "DUF1700-ahelical") %>%
#   str_replace_all("DUF1707-SHOCT-bihelical", "SHOCT-bihelical") %>%
#   str_replace_all("Toast-rack", "Toast_rack") %>%
#   str_replace_all("PADR-HTH", "PadR-HTH") %>%
#   str_replace_all("clgR-HTH", "ClgR-HTH") %>%
#   str_replace_all("ZNR", "ZnR") %>%
#   str_replace_all("Psp-AA", "PspAA") %>%
#   str_replace_all("pspF", "PspF") %>%
#   str_replace_all("Flot", "TM-Flotillin") %>%
#   str_replace_all("phageshock_PspD", "PspD")
#
#
# write_tsv(all, path = "data/rawdata_tsv/all_clean0920.txt")

all <- read_tsv("data/rawdata_tsv/all_clean0920.txt", col_names = T)


domains_of_interest <- c("ZnR", "SIGMA-HTH", "GNTR-HTH", "Beta-Propeller", "DUF2089-HTH", "SHOCT-like",
                         "SHOCT−bihelical",
                         "MacB_PCD", "FTSW_RODA_SPOVE", "PADR−HTH", "ahelical",
                         "REC", "HISKIN", "LiaI−LiaF−TM", "Toast_rack", "PspC", "PspB", "PspF",
                         "Tfu_1009", "PspAA", "Spermine_synth", "TM-Flotillin", "Band_7", "Cest_Tir",
                         "PspA", "Snf7", "Classical-AAA", "TraJ-RHH", "clgR-HTH", "Thioredoxin",
                         "DUF3046","PspN_N", "PspM")
# all$DomArch <- purrr::map(all$DomArch, function(x)str_replace_all(x, pattern = "DUF1707-SHOCT", "DUF1707-SHOCT-bihelical" ))


dom_interest_2 <- c("PspA", "Snf7","Classical-AAA",
                    "PspF","PspB","PspC","ClgR"
                    ,"PspM","Thioredoxin","PspN_N","DUF3046",
                    "LiaI-LiaF-TM", "Toast_rack", "REC", "HISKIN", "ahelical",
                    "SHOCT-bihelical", "SHOCT-like", "Tfu_1009", "PspAA", "Spermine_synth",
                    "TM-Flotillin", "Band_7",
                    "Beta-Propeller", "MacB_PCD", "FTSW_RODA_SPOVE", "Cest_Tir",
                    "SIGMA-HTH", "GNTR-HTH", "DUF2089-HTH", "PadR-HTH",
                    "TraJ-RHH",
                    "ZnR")

domains_of_interest <- dom_interest_2



LineagePlot <- function(prot, domains_of_interest, level = 3)
{
  #' @author Samuel Chen
  #' Generate a lineage plot
  #'
  #' @param prot Data frame containing DomArch and Lineage Columns
  #' @param domains_of_interest Domains to check for the presence of in all the lineages
  #' @param level The max depth of Lineage. ie) i = Kingdom, 2 = Phylum, 3 = class ...



  LevelReduction <- function(lin)
  {
    if(level == 1)
    {
      gt_loc <- str_locate(lin, ">")[[1]]
      if(is.na(gt_loc))
      {
        # No '>' in lineage
        return(lin)
      }
      else
      {
        lin <- substring(lin, first = 0, last = (gt_loc-1))
        return(lin)
      }
    }
    #### Add guard here to protect from out of bounds
    gt_loc <- str_locate_all(lin, ">")[[1]] # [(level-1),][1]
    l <- length(gt_loc)/2
    if( level > l)
    {
      # Not enough '>' in lineage
      return(lin)
    }
    else
    {
      gt_loc <- gt_loc[level,][1] %>% as.numeric()
      lin <- substring(lin, first = 0, last = (gt_loc-1))
      return(lin)
    }
  }

  all_grouped <- data.frame("Query" = character(0), "Lineage" = character(0), "count"= integer())
  for(dom in domains_of_interest)
  {
    domSub <- prot %>% filter(grepl(dom, GenContext, ignore.case = T))

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

  all_grouped <- all_grouped%>%mutate(ReducedLin = unlist(purrr::map(Lineage, LevelReduction)))

  all_grouped_reduced <- all_grouped %>% group_by(Query, ReducedLin) %>% summarize("count"= sum(count)) %>%
    mutate(Kingdom = unlist(purrr::map(ReducedLin,GetKingdom)) )

  lin_counts <- all_grouped_reduced %>% group_by(Kingdom, ReducedLin) %>% summarize("count" =n())

  # grep("bacteria", lin_counts$Kingdom) %>% length()

  bacteria_colors <- rep("#d94f25", length(grep("bacteria", lin_counts$Kingdom)))

  archaea_colors <- rep("#26662d", length(grep("archaea", lin_counts$Kingdom)))

  eukaryota_colors <- rep("#123d99", length(grep("eukaryota", lin_counts$Kingdom)))

  virus_colors <- rep("#538073", length(grep("virus", lin_counts$Kingdom)))
  colors <- append(archaea_colors, bacteria_colors) %>%
    append(eukaryota_colors) %>% append(virus_colors)

  all_grouped_reduced$ReducedLin <- map(all_grouped_reduced$ReducedLin,
                                        function(lin)
                                        {
                                          gt_loc <- str_locate(lin, ">")[,"start"]

                                          if(is.na(gt_loc))
                                          {
                                            # No '>' in lineage
                                            return(lin)
                                          }
                                          else
                                          {
                                            lin <- substring(lin, first = (gt_loc+1), last = 100)
                                            return(lin)
                                          }
                                        }
  ) %>% unlist()

  ordered_lin <- all_grouped_reduced %>% arrange(Kingdom)

  all_grouped_reduced$Query <- factor(x = all_grouped_reduced$Query,
                                      levels = rev(domains_of_interest))
  all_grouped_reduced$ReducedLin <- factor(x = all_grouped_reduced$ReducedLin,
                               levels = unique(ordered_lin$ReducedLin)
  )


  # all_grouped <- all_grouped %>% mutate(Lineage = unlist(purrr::map(Lineage, GetPhylum)))
  # all_grouped_reduced2 <- data.frame("Query" = character(0), "ReducedLin" = character(0), "count"= integer(), "Kingdom" = character(0))
  # for(dom in domains_of_interest)
  # {
  #   domSub <- all_grouped_reduced %>% filter(grepl(dom, Query, ignore.case = T))
  #
  #   all_grouped_reduced2 <- dplyr::union(all_grouped_reduced2, domSub)
  # }
  #
  # all_g




  ggplot(data=all_grouped_reduced,
         aes_string(x="ReducedLin", y="Query")) +
    geom_tile(data=subset(all_grouped_reduced,
                          !is.na(count)),
              aes(fill=count),
              colour="darkred", size=0.3) + #, width=0.7, height=0.7),
    scale_fill_gradient(low="white", high="darkred") +
    # scale_x_discrete(position="top") +
    theme_minimal()+
    theme(axis.title = element_blank() ,axis.text.x=element_text(angle=90,hjust=1,vjust=0.2, color =
                                                                   colors
    ))


}

 LineagePlot(all, domains_of_interest ,2)
ggsave("LineagePlot_2_OldOrder.png",
       width = 10.4,
       height = 5.84,
       dpi = 300
       )

LineagePlot(all, domains_of_interest ,3)
ggsave("LineagePlot_3_OldOrder.png",
       width = 18,
       height = 5.84,
       dpi = 300
)

# LineagePlot(all, domains_of_interest ,10)
# ggsave("LineagePlot_FullLin.png",
#        width = 12.3,
#        height = 5.84,
#        dpi = 300
# )
