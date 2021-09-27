## To convert IPRScan files to a gggenes viz!
## Janani Ravi, Lauren Sosinski, Samuel Chen
## Created: Apr 9, 2020

library(here)
library(tidyverse)
library(data.table)
library(gggenes)
library(ggplot2)
#source("../the-approach/R/pre-msa-tree.R") # for "to_titlecase()"
source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")

#################################
## Modified gggenes::theme_genes
#################################
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
    strip.background = ggplot2::element_blank()
    #strip.text = ggplot2::element_blank()
  )
}
CPCOLS <- c('#AFEEEE', '#DDA0DD', '#EE2C2C', '#CDBE70', '#B0B099',
             '#8B2323', '#EE7600', '#EEC900', 'chartreuse3', '#0000FF',
             '#FFD900', '#32CD32', 'maroon4', 'cornflowerblue', 'darkslateblue',
             '#AB82FF', '#CD6889', '#FFA07A', '#FFFF00', '#228B22',
             '#FFFFE0', '#FFEC8B', 'peru', '#668B8B', 'honeydew',
             '#A020F0', 'grey', '#8B4513', '#191970', '#00FF7F',
             'lemonchiffon','#66CDAA', '#5F9EA0', '#A2CD5A', '#556B2F',
             '#EEAEEE', 'thistle4', '#473C8B', '#FFB6C1', '#8B1C62',
             '#FFE4B5', 'black', '#FF7F50', '#FFB90F', '#FF69B4', '#836FFF',
             '#757575','#CD3333', '#EE7600', '#CDAD00', '#556B2F', '#7AC5CD')
##################################
## Get Top N AccNum by Lin+DomArch
##################################
# Group by lineage + DA then take top 20
find_top_acc = function(infile_full,
                        DA_col = "DomArch.Pfam", ## @SAM, you could pick by the Analysis w/ max rows!
                        lin_col = "Lineage",
                        n = 20)
{
  lin_sym = sym(lin_col)
  DA_sym = sym(DA_col)

  cln = fread(infile_full, sep ="\t", fill = T)

  ## Group by Lineage, DomArch and reverse sort by group counts
  grouped = cln %>%
    group_by({{lin_sym}}, {{DA_sym}}) %>%
    summarise(count = n()) %>%
    arrange(-count) %>%
    filter(!is.na({{lin_sym}}) & !is.na({{DA_sym}}))

  top_acc = character(n)
  for(r in 1:min(nrow(grouped), n))
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


#############################################
## IPR + FULL files --> DomArch Visualization
#############################################
ipr2viz <- function(infile_ipr=NULL, infile_full=NULL,
                    analysis=c("Pfam", "Phobius", "TMHMM", "Gene3D"),
                    group_by = "Query", #"Analysis"
                    topn = 20, name = "Name", text_size = 10)
{
  ## Read IPR file
  ipr_out <- read_tsv(infile_ipr, col_names=T)

  ## To filter by Analysis
  analysis = paste(analysis, collapse = "|")

  ## @SAM: This can't be set in stone since the analysis may change!
  ## Getting top n accession numbers using find_top_acc()
  top_acc <- find_top_acc(infile_full=infile_full,
                          DA_col = "DomArch.Pfam",
                          ## @SAM, you could pick by the Analysis w/ max rows!
                          lin_col = "Lineage",
                          n = topn)

  # Filter by Top Accessions per Accession per DomArch and Lineage
  ipr_out <- subset(ipr_out,
                    ipr_out$AccNum %in% top_acc)

  ## Need to fix this eventually based on the 'real' gene orientation! :)
  ipr_out$Strand <- rep("forward", nrow(ipr_out))

  ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out,
                        grepl(pattern=analysis, x=Analysis))

  # dynamic analysis labeller
  analyses <- ipr_out_sub %>%
    select(Analysis) %>%
    distinct()
  analysis_labeler <- analyses %>%
    pivot_wider(names_from = Analysis, values_from = Analysis)
  print(analysis_labeler)

  queryrows <- which(is.na(ipr_out_sub$AccNum))

  lookup_tbl_path = "/data/research/jravilab/common_data/cln_lookup_tbl.tsv"
  lookup_tbl = read_tsv(lookup_tbl_path, col_names = T)

  print(colnames(ipr_out_sub))
  print(colnames(lookup_tbl))
  lookup_tbl = lookup_tbl %>% select(-ShortName) # Already has ShortName -- Just needs SignDesc
  # ipr_out_sub = ipr_out_sub %>% select(-ShortName)
  ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by = "DB.ID")

  print(colnames(ipr_out_sub))
  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "Analysis")
  {
    ggplot(ipr_out_sub,
           aes_string(xmin = "StartLoc", xmax = "StopLoc",
                      y = name, fill = "SignDesc", label="ShortName")) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      #geom_blank(data = dummies) +
      facet_wrap(~ Analysis, strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      #, ncol = 1 + #scales = "free",
      scale_color_manual(values = CPCOLS) +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

  else if(group_by == "Query"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis,  #y = AccNum
               fill = SignDesc, label = ShortName)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      facet_wrap(as.formula(paste("~", name)), strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_color_manual(values = CPCOLS)  +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

}

ipr2viz_web <- function(infile_ipr,
                        accessions,
                        analysis=c("Pfam", "Phobius","TMHMM","Gene3D"),
                        group_by = "Query", name = "Name",
                        text_size = 10)
{
  ## To filter by Analysis
  analysis = paste(analysis, collapse = "|")

  ## @SAM, colnames, merges, everything neeeds to be done now based on the
  ## combined lookup table from "common_data"
  lookup_tbl_path = "/data/research/jravilab/common_data/cln_lookup_tbl.tsv"
  lookup_tbl = read_tsv(lookup_tbl_path, col_names = T)

  ## Read IPR file and subset by Accessions
  ipr_out <- read_tsv(infile_ipr, col_names=T)
  #ipr_out <- subset(ipr_out, ipr_out$AccNum %in% accessions)

  ## Need to fix eventually based on 'real' gene orientation!
  ipr_out$Strand <- rep("forward", nrow(ipr_out))

  ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out,
                        grepl(pattern=analysis, x=Analysis))

  # dynamic analysis labeller
   analyses <- ipr_out_sub %>%
    select(Analysis) %>%
    distinct()
  analysis_labeler <- analyses %>%
    pivot_wider(names_from = Analysis, values_from = Analysis)
  # analysis_labeler[1,] = colnames(analysis_labeler)

  print(analysis_labeler)

  queryrows <- which(is.na(ipr_out_sub$AccNum))
  lookup_tbl = lookup_tbl %>% select(-ShortName)
  ## @SAM, make sure the following two work with the Lookup Tables!!
  ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by = "DB.ID")

  ## PLOTTING
  ## domains as separate arrows
  
  if(group_by == "Analysis")
  {
    ggplot(ipr_out_sub,
           aes_string(xmin = "StartLoc", xmax = "StopLoc",
                      y = name, fill = "SignDesc", label="ShortName")) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      #geom_blank(data = dummies) +
      facet_wrap(~ Analysis, strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      #, ncol = 1 + #scales = "free",
      scale_color_manual(values = CPCOLS) +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

  else if(group_by == "Query"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis,  #y = AccNum
               fill = SignDesc, label=ShortName)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      facet_wrap(as.formula(paste("~", name)),
                 strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_color_manual(values = CPCOLS) +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=3))
  }

}
