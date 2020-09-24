## To convery IPRScan files to a gggenes viz!
## Janani Ravi, Lauren Sosinski, Samuel Chen
## Created: Apr 9, 2020

library(here)
library(tidyverse)
library(gggenes)
library(ggplot2)
#source("../the-approach/R/pre-msa-tree.R") # for "to_titlecase()"

## theme_genes2 adapted from theme_genes (w/o strip.text())
## https://github.com/wilkox/gggenes/blob/master/R/theme_genes.R
theme_genes2 <- function() {
  ggplot2::theme_grey() + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey", size = 1),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
    axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
    #strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()
  )
}

ipr2domarch <- function(infile_ipr=here("evolvr/TestData/Iprscan/sausa300_0204.iprscan.tsv"),
                        PfamClans_path = "evolvr/TestData/Pfam-A.clans.txt",
                        analysis=c("Pfam", "Phobius","TMHMM","SUPERFAMILY"),
                        group_by = "DB",
                        topn = 20)
{

  analysis = paste(analysis, collapse = "|")

  colm <- c("SignAcc", "ClID", "Short_Name", "Name", "Description")

  PfamClans = read_tsv(PfamClans_path, col_names = colm)

  ## Read in Files and Assign ColNames
  ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
                "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
                "Status", "RunDate",
                "IPRAcc", "IPRDesc")#, "GOAnn", "PathAnn") ## OPTIONAL
  ipr_out <- read_tsv(infile_ipr, col_names=ipr_cols, skip = 1) #%>%
  #left_join(protmap, by="AccNum")

  ## getting top n accession numbers
  top_acc <- unique(ipr_out$AccNum) %>%
    head(topn)
  ipr_out <- subset(ipr_out, ipr_out$AccNum %in% top_acc)


  #colnames(ipr_out) <- ipr_cols
  ipr_out$Strand <- rep("forward", nrow(ipr_out))


  #########ipr_out <- ipr_out %>% arrange(Leaf, StartLoc, StopLoc)
  ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out, grepl(pattern=analysis, x=Analysis))


  # dynamic analysis labeller
  analyses <- ipr_out_sub %>% select(Analysis) %>%
    distinct()
  analysis_labeler <- analyses %>% pivot_wider(names_from = Analysis, values_from = Analysis)

  queryrows <- which(is.na(ipr_out_sub$AccNum))
  # ipr_out_sub$Leaf[queryrows] = ipr_out_sub$AccNum[queryrows]
  #ipr_out_sub$ProtName[queryrows] #"sausa300_0204"

  ipr_out_sub <- left_join(ipr_out_sub, PfamClans, by = "SignAcc")

  signacc <- which(is.na(ipr_out_sub$Name))
  ipr_out_sub$Name[signacc] = ipr_out_sub$SignAcc[signacc]

  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "DB")
  {
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = AccNum, fill = SignDesc, label=Name)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +
      #geom_blank(data = dummies) +
      facet_wrap(~ Analysis, strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      #, ncol = 1 + #scales = "free",
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom", legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin()) +
      ylab("")
  }

  else if(group_by == "protein"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis  #y = AccNum
               , fill = SignDesc, label=Name)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +

      facet_wrap(~ AccNum, strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom", legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin()) +
      ylab("")
  }

}
