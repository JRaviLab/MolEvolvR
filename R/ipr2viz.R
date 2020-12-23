## To convert IPRScan files to a gggenes viz!
## Janani Ravi, Lauren Sosinski, Samuel Chen
## Created: Apr 9, 2020

library(here)
library(tidyverse)
library(data.table)
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

# Group by lineage + DA then take top 20
top_acc = function(cln_file, DA_col = "DomArch.Pfam", lin_col = "Lineage", n = 20)
{
  lin_sym = sym(lin_col)
  DA_sym = sym(DA_col)

  cln = fread(cln_file, sep ="\t", fill = T)

  grouped = cln %>% group_by({{lin_sym}}, {{DA_sym}}) %>% summarise(count = n()) %>%
    arrange(-count) %>% filter(!is.na({{lin_sym}}) & !is.na({{DA_sym}}))

  top_acc = character(n)
  for(r in 1:min(nrow(grouped),n))
  {
    l = (grouped %>% pull({{lin_sym}}))[r]
    DA = (grouped %>% pull({{DA_sym}}))[r]

    filt = cln %>% filter({{lin_sym}} == l & {{DA_sym}} == DA)

    top = filt[which(filt$PcPositive == max(filt$PcPositive))[1] , ]

    top_acc[r] = top$AccNum
  }
  top_acc = top_acc[which(top_acc != "")]
  return(top_acc)
}






ipr2viz <- function(infile_ipr=here("../molevolvr_out/phage1_out/WP_001901328.1_Vibrio_cholerae.iprscan.tsv"),
                    cln_file = "../molevolvr_out/phage1_out/cln_combined.tsv",
                    PfamClans_path = "TestData/Pfam-A.clans.txt",
                    analysis=c("Pfam", "Phobius","TMHMM","SUPERFAMILY"),
                    group_by = "DB",
                    topn = 20, name = "AccNum", text_size = 10)
{

  analysis = paste(analysis, collapse = "|")

  colm <- c("SignAcc", "ClID", "Short_Name", "PFC_Name", "Description")

  PfamClans = read_tsv(PfamClans_path, col_names = colm)

  ## Read in Files and Assign ColNames
  ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
                "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
                "Status", "RunDate",
                "IPRAcc", "IPRDesc")#, "GOAnn", "PathAnn") ## OPTIONAL
  ipr_out <- read_tsv(infile_ipr, col_names=ipr_cols, skip = 1) #%>%
  #left_join(protmap, by="AccNum")

  ## getting top n accession numbers
  top_acc <- top_acc(cln_file, DA_col = "DomArch.Pfam", lin_col = "Lineage",
                     n = topn)
  # Add the lineages from the cln file
  # lin_cols = c("AccNum", "Lineage", "Species")
  ipr_out = merge(ipr_out, cln_file[, c("AccNum","Lineage", "Species")], by = "AccNum", all.x = T)


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
  print(analysis_labeler)

  queryrows <- which(is.na(ipr_out_sub$AccNum))

  ipr_out_sub <- left_join(ipr_out_sub, PfamClans, by = "SignAcc")

  signacc <- which(is.na(ipr_out_sub$PFC_Name))
  ipr_out_sub$PFC_Name[signacc] = ipr_out_sub$SignAcc[signacc]

  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "DB")
  {
    ggplot(ipr_out_sub,
           aes_string(xmin = "StartLoc", xmax = "StopLoc",
               y = name, fill = "SignDesc", label="PFC_Name")) +
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
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

  else if(group_by == "protein"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis  #y = AccNum
               , fill = SignDesc, label=PFC_Name)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +

      facet_wrap(as.formula(paste("~", name)), strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom", legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

}

ipr2viz_noblast <- function(infile_ipr,
                    accessions,
                    PfamClans_path = "TestData/Pfam-A.clans.txt",
                    analysis=c("Pfam", "Phobius","TMHMM","SUPERFAMILY"),
                    group_by = "DB", name = "AccNum", text_size = 10)
{

  analysis = paste(analysis, collapse = "|")

  colm <- c("SignAcc", "ClID", "Short_Name", "PFC_Name", "Description")

  PfamClans = read_tsv(PfamClans_path, col_names = colm)

  ## Read in Files and Assign ColNames
  # ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
  #               "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
  #               "Status", "RunDate",
  #               "IPRAcc", "IPRDesc")#, "GOAnn", "PathAnn") ## OPTIONAL
  ipr_out <- fread(infile_ipr, fill = T, sep = "\t")
  #left_join(protmap, by="AccNum")

  ipr_out <- subset(ipr_out, ipr_out$AccNum %in% accessions)


  #colnames(ipr_out) <- ipr_cols
  ipr_out$Strand <- rep("forward", nrow(ipr_out))


  #########ipr_out <- ipr_out %>% arrange(Leaf, StartLoc, StopLoc)
  ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out, grepl(pattern=analysis, x=Analysis))

  # dynamic analysis labeller
  analyses <- ipr_out_sub %>% select(Analysis) %>%
    distinct()
  analysis_labeler <- analyses %>% mutate(id = 1) %>% pivot_wider(names_from = Analysis,
                                                             values_from = Analysis) %>%
    select(-id)
  # analysis_labeler[1,] = colnames(analysis_labeler)

  print(analysis_labeler)

  queryrows <- which(is.na(ipr_out_sub$AccNum))

  ipr_out_sub = dplyr::rename(ipr_out_sub,  "SignAcc" = "DB.ID")
  ipr_out_sub <- left_join(ipr_out_sub, PfamClans, by = "SignAcc")

  signacc <- which(is.na(ipr_out_sub$PFC_Name))
  ipr_out_sub$PFC_Name[signacc] = ipr_out_sub$SignAcc[signacc]

  name_sym = sym(name)
  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "DB")
  {
    ggplot(ipr_out_sub,
           aes_string(xmin = "StartLoc", xmax = "StopLoc",
               y = name, fill = "SignDesc", label="PFC_Name")) +
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
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

  else if(group_by == "protein"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis  #y = AccNum
               , fill = SignDesc, label=PFC_Name)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +

      facet_wrap(as.formula(paste("~", name)), strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom", legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"), legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

}
