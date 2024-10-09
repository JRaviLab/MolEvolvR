library(tidyverse)
library(aplot)
library(ape)
library(tidytree)
library(seqinr)
library(msa)
library(htmltools)
library(data.table)

# figures 2 and 4 from DciA paper, https://journals.asm.org/doi/epub/10.1128/jb.00163-22

df_shared <- function(job_dir) {
  df <- read_tsv(file.path(job_dir, "./query_data/query_data.ipr_domarch.tsv")) # /data/scratch/jburke/dcia/all_seqs_full_with_mobi_ipr_distinct_2.tsv

  df$Name <- str_replace_all(df$Name, "\\[", "")
  df$QueryName <- df$Name

  df <- df %>% filter(Lineage != "") %>% drop_na(Lineage)
  df <- df %>% mutate(DomArch.Pfam = gsub("_[0-9]{1,2}", "", DomArch.Pfam)) %>%
    mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "TPR\\+.*", paste0("TPR(", str_count(DomArch.Pfam, "TPR"), ")"))) %>%
    mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "DciA\\+DciA\\+DciA\\+DciA", "DciA(4)"))

  return(df)
}

dcia_figure2 <- function(job_dir) {
  ## Figure 2
  # Fig 2a length box plot of DciA only bac
  
  df <- df_shared(job_dir)

  df_fig2 <- df %>%
    filter(DomArch.Pfam == "DciA") %>%
    filter(grepl("Bacteria>", Lineage)) %>%
    filter(!grepl("Candidatus", Lineage)) %>%
    filter(!grepl("candidate", Lineage)) %>%
    distinct(AccNum, .keep_all = TRUE) %>%
    drop_na(pfam_start)

  df_fig2$domainLength <- df_fig2$pfam_stop - df_fig2$pfam_start
  df_fig2$fromStart <- df_fig2$pfam_start
  df_fig2$fromEnd <- df_fig2$seq_length - df_fig2$pfam_stop
  df_fig2$seqLength <- df_fig2$seq_length
  df_fig2 <- df_fig2 %>%
    pivot_longer(c(seqLength, fromStart, fromEnd), names_to = "Type", values_to = "length")
  df_fig2$Type <- factor(df_fig2$Type, levels = c("seqLength", "fromStart", "fromEnd"))
  df_fig2_stats <- df_fig2 %>%
    group_by(Type) %>%
    mutate(med = median(length), twent = quantile(length, c(0.25)), sev = quantile(length, c(0.75)))

  fig_2a <- ggplot(df_fig2, aes(x=Lineage, length)) +
    geom_boxplot(outlier.alpha = 0.5, outlier.color = "black") +
    labs(y = "Length", x = "Lineage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.background = element_rect("white", size = 0),
          text = element_text(size = 18), panel.border = element_blank()) +
    facet_grid(Type ~ . , scales = "free") +
    geom_hline(data = df_fig2_stats, aes(yintercept = med), color = "blue") +
    geom_hline(data = df_fig2_stats, aes(yintercept = twent), color = "blue", alpha = 0.5) +
    geom_hline(data = df_fig2_stats, aes(yintercept = df_fig2_stats$sev), color = "blue", alpha = 0.5)

  return(fig_2a)
}

dcia_figure4 <- function(job_dir, save_figs=FALSE) {
  df <- df_shared(job_dir)

  df_bact <- df %>% filter(grepl("Bacteria>", Lineage))

  df_bact <- df_bact %>%
    filter(grepl("Bacteria", Lineage)) %>%
    drop_na(AccNum) %>%
    filter(!grepl("candidate", Lineage)) %>%
    filter(!grepl("Candidatus", Lineage)) %>%
    filter(!grepl("candidate", Species)) %>%
    filter(!grepl("Candidat", Species)) %>%
    filter(!grepl("uncultur", Species))

  df_bact$DomArch.Pfam <- str_replace_all(df_bact$DomArch.Pfam , "Hydrolase\\+DciA", "HAD_2+DciA")

  ## FIGURE 4
  # Fig 4a stacked bar plot of lineage
  df_fig4 <- df_bact %>% filter(DomArch.Pfam == "DciA")%>% drop_na(seq_length)
  df_fig4$fromEnd <- df_fig4$seq_length - df_fig4$pfam_stop
  df_fig4$fromStart <- df_fig4$pfam_start
  df_fig4$group <- ""
  df_fig4$group[which(df_fig4$fromStart >= 14 & df_fig4$fromEnd < 14)] <- 1
  df_fig4$group[which(df_fig4$fromEnd >= 14 &  df_fig4$fromStart < 14)] <- 2
  df_fig4$group[which(df_fig4$fromEnd >= 14 &  df_fig4$fromStart >= 14)] <- 3
  df_fig4$group[which(df_fig4$group == "")] <- 4
  write_tsv(df_fig4, "dcia_figs/fig4/cln_with_groups.tsv")

  stacked_fig4a <- stacked_lin_plot(df_fig4, column = "group", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.3),
                              legend.cols = 2)

  stacked_fig4a_nolegend <- stacked_lin_plot(df_fig4, column = "group", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.3),
                              legend.cols = 2, legend = FALSE)

  # Fig 4b
  df_fig4_b <- df_fig4 %>% filter(Lineage == "Bacteria>Proteobacteria" | Lineage == "Bacteria>Actinobacteria" | Lineage == "Bacteria>Bacteroidetes")
  stacked_fig4b <- stacked_lin_plot(df_fig4_b, column = "group",
                              Lineage_col = "Lineage_short", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.4),
                              legend.cols = 2, reduce_lineage = FALSE)
                              Lineage_col = "Lineage_short", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.4),
                              legend.cols = 2, reduce_lineage = FALSE, legend = FALSE)

  # Fig 4c
  df_fig4_c <- df_fig4 %>%
    filter(Lineage != "Bacteria>Proteobacteria" & Lineage != "Bacteria>Actinobacteria" & Lineage != "Bacteria>Bacteroidetes")
  stacked_fig4c <- stacked_lin_plot(df_fig4_c, column = "group", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.15),
                              legend.cols = 4)
  stacked_fig4c_nolegend <- stacked_lin_plot(df_fig4_b, column = "group", cutoff = 100,
                              xlabel = "Group", legend.position = c(0.7,0.15),
                              legend.cols = 4, legend = FALSE)

  if (save_figs) {
    ggsave("dcia_figs/fig4/4a.png", stacked_fig4a, dpi = 400, device = "png", height = 15, width = 15)
    ggsave("dcia_figs/fig4/4_no_legend.png", stacked_fig4a_nolegend, dpi = 400, device = "png", height = 15, width = 15)
    
    ggsave("dcia_figs/fig4/4b.png", stacked_fig4b, dpi = 400, device = "png", height = 15, width = 15) stacked_fig4b_nolegend <- stacked_lin_plot(df_fig4_b, column = "group",
    ggsave("dcia_figs/fig4/4b_no_legend.png", stacked_fig4b_nolegend, dpi = 400, device = "png", height = 15, width = 15)
    
    ggsave("dcia_figs/fig4/4c.png", stacked_fig4c, dpi = 400, device = "png", height = 15, width = 15)
    ggsave("dcia_figs/fig4/4c_no_legend.png", stacked_fig4c_nolegend, dpi = 400, device = "png", height = 15, width = 15)
  }

  return(
    list(
      stacked_fig4a = stacked_fig4a,
      stacked_fig4a_nolegend = stacked_fig4a_nolegend,

      stacked_fig4b = stacked_fig4b,
      stacked_fig4b_nolegend = stacked_fig4b_nolegend,
      
      stacked_fig4c = stacked_fig4c,
      stacked_fig4c_nolegend = stacked_fig4c_nolegend
    )
  )
}
