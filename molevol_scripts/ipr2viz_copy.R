## To convert IPRScan files to a gggenes viz!
## Janani Ravi, Lauren Sosinski, Samuel Chen
## Created: Apr 9, 2020

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gggenes))
suppressPackageStartupMessages(library(ggplot2))
#source("../the-approach/R/pre-msa-tree.R") # for "to_titlecase()"
source("R/colnames_molevol.R")
source("R/combine_files.R")
source("R/combine_analysis.R")

#################################
## Modified gggenes::theme_genes
#################################
## theme_genes2 adapted from theme_genes (w/o strip.text())
## https://github.com/wilkox/gggenes/blob/master/R/theme_genes.R
theme_genes2 <- function() {
  ggplot2::theme_grey() + ggplot2::theme(
    panel.background=ggplot2::element_blank(),
    panel.grid.major.y=ggplot2::element_line(colour="grey", size=1),
    panel.grid.minor.y=ggplot2::element_blank(),
    panel.grid.minor.x=ggplot2::element_blank(),
    panel.grid.major.x=ggplot2::element_blank(),
    axis.ticks.y=ggplot2::element_blank(),
    axis.line.x=ggplot2::element_line(colour="grey20", size=0.5),
    axis.ticks.x=ggplot2::element_line(colour="grey20", size=0.5),
    strip.background=ggplot2::element_blank()
    #strip.text=ggplot2::element_blank()
  )
}

##################################
## Get Top N AccNum by Lin+DomArch
##################################
# Group by lineage + DA then take top 20
find_top_acc <- function(infile_full,
                         DA_col="DomArch.Pfam", ## @SAM, you could pick by the Analysis w/ max rows!
                         lin_col="Lineage",
                         n=20)
{
  lin_sym <- sym(lin_col)
  DA_sym <- sym(DA_col)

  cln <- fread(infile_full, sep ="\t", fill=T)

  ## Group by Lineage, DomArch and reverse sort by group counts
  grouped <- cln %>%
    group_by({{lin_sym}}, {{DA_sym}}) %>%
    summarise(count=n()) %>%
    arrange(-count) %>%
    filter(!is.na({{lin_sym}}) & !is.na({{DA_sym}}))

  top_acc <- character(n)
  for(r in 1:min(nrow(grouped), n))
  {
    l <- (grouped %>% pull({{lin_sym}}))[r]
    DA <- (grouped %>% pull({{DA_sym}}))[r]

    filt <- cln %>% filter({{lin_sym}} == l & {{DA_sym}} == DA)

    top <- filt[which(filt$PcPositive == max(filt$PcPositive))[1] , ]

    top_acc[r] <- top$AccNum
  }
  top_acc <- top_acc[which(top_acc != "")]
  return(top_acc)
}

#############################################
## IPR + FULL files --> DomArch Visualization
#############################################
ipr2viz <- function(infile_ipr=NULL, infile_full=NULL,
                    analysis=c("Pfam", "Gene3D", "Phobius"),
                    group_by="Analysis",
                    topn=20, name="Name", text_size=12)
{

  ## Populating function ARGs temporarily
  group_by='Analysis' #; group_by='Query'
  # analysis="ProSiteProfiles"; analysis="Pfam"; analysis="PANTHER"
  analysis <- c("Pfam", "Gene3D", "PANTHER",
                "Phobius", "TMHMM", #"SignalP_GRAM_POSITIVE",
                "SUPERFAMILY",
                "ProSiteProfiles", "MobiDBLite")
  topn=F; name="Name"; text_size=14

  ##!! @LAUREN's special case for Staph
  # prefix = 'Saureus_gisC_ABD20640.1'

  ## EXAMPLE DATASETS for plotting ##
  ## SLPS
  # inpath <- "../molevol_data/project_data/slps/da_analysis_20210116/"
  # infile_ipr <- paste0(inpath, 'ipr_combined.tsv', collapse="")
  ## PHAGE DEFENSE
  # inpath <- "../molevol_data/project_data/phage_defense/full_analysis_20210108/"
  # infile_ipr <- paste0(inpath, 'dcdv_quick_out/dcdv.iprscan_cln.tsv',
  #                      collapse="")
  # infile_full <- '../full_analysis_20210108/WP_001901328_full/WP_001901328.full_analysis.tsv'
  ## STAPH INPUTS
  # infile_ipr <- '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.iprscan_cln.tsv'
  # infile_full <- '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.full_analysis.tsv'

  ## Read IPR file
  ipr_out <- read_tsv(infile_ipr, col_names=T)

  ##!! @LAUREN's special case for Staph
  # query <- ipr_out %>%
  #   filter(AccNum == 'ABD20640.1')
  # ipr_out <- ipr_out %>%
  #   filter(grepl(pattern="BFirmi|BActi", Name))

  ## To filter by Analysis
  analysis <- paste(analysis, collapse="|")

  ## If filtered AccNum need to be visualized
  if((topn > 0) && (topn != F) && (topn != FALSE)){
    ## @SAM: This can't be set in stone since the analysis may change!
    ## Getting top n accession numbers using find_top_acc()
    top_acc <- find_top_acc(infile_full=infile_full,
                            DA_col="DomArch.Pfam",
                            ## @SAM, you could pick by the Analysis w/ max rows!
                            lin_col="Lineage",
                            n=topn)

    # Filter by Top Accessions per Accession per DomArch and Lineage
    ipr_out <- subset(ipr_out,
                      ipr_out$AccNum %in% top_acc)
  }

  ## Need to fix this eventually based on the 'real' gene orientation! :)
  ipr_out$Strand <- rep("forward", nrow(ipr_out))
  ##!! @LAUREN's special case for Staph
  # query$Strand <- rep("forward", nrow(query))

  ## Order domains from Start -> End within a protein
  ipr_out <- ipr_out %>%
    arrange(AccNum, StartLoc, StopLoc) %>%
    group_by(AccNum, Analysis, StartLoc, StopLoc) %>%
    slice_head(n=1)

  # Predominant analysis for this IPRSCAN run
  table(ipr_out$Analysis) %>%
    sort(decreasing=T)

  ##!! @LAUREN's special case for Staph
  # ipr_out_sub <- query %>%
    # bind_rows(ipr_out)

  # Subset by 'pre-selected analysis'
  ipr_out_sub <- ipr_out %>%
    filter(grepl(pattern=analysis, x=Analysis))

  # dynamic analysis labeller
  analyses <- ipr_out_sub %>%
    select(Analysis) %>%
    distinct()

  analysis_labeler <- analyses %>%
    pivot_wider(names_from=Analysis, values_from=Analysis)
  print(analysis_labeler)

  # Check for missing AccNum
  # queryrows <- which(is.na(ipr_out_sub$AccNum))
  # cat("Missing AccNum indices: ", queryrows)

  ## Custom color palette
  CPCOLS <- c('#AFEEEE', '#DDA0DD', '#EE2C2C', '#CDBE70', '#B0B099',
              '#8B2323', '#EE7600', '#EEC900', 'chartreuse3', '#0000FF',
              '#FFD900', '#32CD32', 'maroon4', 'cornflowerblue', #'darkslateblue',
              '#AB82FF', '#CD6889', '#FFA07A', '#FFFF00', '#228B22',
              '#FFFFE0', '#FFEC8B', 'peru', '#668B8B', 'honeydew',
              '#A020F0', 'grey', '#8B4513', '#191970', '#00FF7F',
              'lemonchiffon','#66CDAA', '#5F9EA0', '#A2CD5A', '#556B2F',
              '#EEAEEE', 'thistle4', '#473C8B', '#FFB6C1', '#8B1C62',
              '#FFE4B5', 'black', '#FF7F50', '#FFB90F', '#FF69B4', '#836FFF',
              '#757575','#CD3333', '#EE7600', '#CDAD00', '#556B2F', '#7AC5CD')

  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "Analysis"){
    ggplot(ipr_out_sub,
           aes_string(xmin = "StartLoc", xmax = "StopLoc",
                      y = name, fill = "SignDesc")) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"),
                      arrow_body_height = unit(3, "mm")) +
      #geom_gene_label(align = "left", min.size=6, grow=T, reflow=T) +
      #geom_blank(data = dummies) +
      facet_wrap(~ Analysis, strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      #, ncol = 1 + #scales = "free",
      # scale_color_manual(values = CPCOLS) +
      scale_fill_brewer(palette="Set2") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            text = element_text(size = 13)) +
      ylab("")+
      guides(fill=guide_legend(nrow=4))
  }

  else if(group_by == "Query"){
    ggplot(ipr_out_sub,
           aes(xmin = StartLoc, xmax = StopLoc,
               y = Analysis,  #y = AccNum
               fill = SignDesc, label = Label)) +
      geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm")) +
      geom_gene_label(align = "left") +

      facet_wrap(as.formula(paste("~", name)), strip.position = "top", ncol = 3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_fill_manual(values = CPCOLS) +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box = "horizontal",
            legend.key.size = unit(0.02, "npc"),
            legend.box.margin = margin(),
            legend.text = element_text(size = 7),
            legend.title = element_blank(),
            text = element_text(size = text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=4))
  }
}

ipr2viz_web <- function(infile_ipr,
                        accessions,
                        analysis=c("Pfam", "Phobius","Gene3D"),
                        group_by="Analysis", name="AccNum",
                        text_size=10)
{
  ## To filter by Analysis
  analysis <- paste(analysis, collapse="|")

  ## @SAM, colnames, merges, everything neeeds to be done now based on the
  ## combined lookup table from "common_data"
  lookup_tbl_path <- "/data/research/jravilab/common_data/lookup_tbl.tsv"
  lookup_tbl <- read_tsv(lookup_tbl_path, col_names=T)

  ## Read IPR file and subset by Accessions
  ipr_out <- read_tsv(infile_ipr, col_names=T)
  ipr_out <- subset(ipr_out, ipr_out$AccNum %in% accessions)

  ## Need to fix eventually based on 'real' gene orientation!
  ipr_out$Strand <- rep("forward", nrow(ipr_out))

  ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out,
                        grepl(pattern=analysis, x=Analysis))

  # dynamic analysis labeller
  analyses <- ipr_out_sub %>% select(Analysis) %>%
    distinct()
  analysis_labeler <- analyses %>% mutate(id=1) %>%
    pivot_wider(names_from=Analysis,
                values_from=Analysis) %>%
    select(-id)
  # analysis_labeler[1,]=colnames(analysis_labeler)

  print(analysis_labeler)

  queryrows <- which(is.na(ipr_out_sub$AccNum))

  ## @SAM, make sure the following two work with the Lookup Tables!!
  ## These two rows aren't needed, lookup table being used doesn't
  ## have a SignAcc column, nor do iprscan results anymore -- all
  ## been changed to DB.ID

  # ipr_out_sub=dplyr::rename(ipr_out_sub,  "SignAcc"="DB.ID")
  # ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by="SignAcc")

  signacc <- which(is.na(ipr_out_sub$ShortName))
  ipr_out_sub$ShortName[signacc] <- ipr_out_sub$DB.ID[signacc]

  name_sym <- sym(name)

  ## PLOTTING
  ## domains as separate arrows
  if(group_by == "Analysis")
  {
    ggplot(ipr_out_sub,
           aes_string(xmin="StartLoc", xmax="StopLoc",
                      y=name, fill="SignDesc", label="ShortName")) +
      geom_gene_arrow(arrowhead_height=unit(3, "mm"),
                      arrowhead_width=unit(1, "mm")) +
      geom_gene_label(align="left") +
      #geom_blank(data=dummies) +
      facet_wrap(~ Analysis, strip.position="top", ncol=3,
                 labeller=as_labeller(analysis_labeler)) +
      #, ncol=1 + #scales="free",
      scale_fill_brewer(palette="Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box="horizontal",
            legend.key.size=unit(0.02, "npc"),
            legend.box.margin=margin(),
            text=element_text(size=text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=4))
  }

  else if(group_by == "Query"){
    ggplot(ipr_out_sub,
           aes(xmin=StartLoc, xmax=StopLoc,
               y=Analysis,  #y=AccNum
               fill=SignDesc, label=ShortName)) +
      geom_gene_arrow(arrowhead_height=unit(3, "mm"),
                      arrowhead_width=unit(1, "mm")) +
      geom_gene_label(align="left") +

      facet_wrap(as.formula(paste("~", name)),
                 strip.position="top", ncol=3,
                 labeller=as_labeller(analysis_labeler)) +
      scale_fill_brewer(palette="Set3") +
      theme_minimal() + theme_genes2() +
      theme(legend.position="bottom",
            legend.box="horizontal",
            legend.key.size=unit(0.02, "npc"),
            legend.box.margin=margin(),
            text=element_text(size=text_size)) +
      ylab("")+
      guides(fill=guide_legend(nrow=4))
  }

}
