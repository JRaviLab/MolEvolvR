## To convery IPRScan files to a gggenes viz!
## Janani Ravi
## Created: Apr 9, 2020

library(here)
library(tidyverse)
library(gggenes)
source("../the-approach/R/pre-msa-tree.R") # for "to_titlecase()"

ipr2domarch <- function(infile_ipr=here("data/saureus/sausa300_0204.1stprotclust89.iprscan.tsv"),
                        infile_blast=here("data/saureus/sausa300_0204.nr.1e5.txt"),
                        analysis="Pfam|Phobius|TMHMM")
{
  #infile=here("data/test_data/sample.iprscan5.tsv"); sppname="Mtb"; analysis="Pfam"
  #infile=here("data/banthracis/ban_sap.iprscan5.tsv"); sppname="B. anthracis"; analysis="Pfam"
  #infile_ipr=here("data/test_data/sample.iprscan5.tsv"); infile_blast=here("data/test_data/sample.nr.1e5.txt"); analysis="Pfam|Phobius|TMHMM"


  ## Read in Files and Assign ColNames
  ipr_cols <- c("AccNum", "Seq_MD5_digest", "SeqLen", "Analysis",
                "SignAcc", "SignDesc", "StartLoc", "StopLoc", "Score",
                "Status", "RunDate",
                "IPRAcc", "IPRDesc", "GOAnn", "PathAnn") ## OPTIONAL
  ipr_out <- read_tsv(infile_ipr, col_names=ipr_cols)
  # based on most recent June 2020, deltablast outfmt
  blast_cols <- c("qacc", "sacc", "sseqid", "sallseqid", "stitle",
                  "sscinames", "staxids", "sskingdoms",
                  "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send",
                  "evalue", "bitscore", "positive", "ppos",
                  "slen", "sgi", "sallgi", "qcovs", "qcovhsp")
  blast_out <- read_tsv(infile_blast,
                        col_names=blast_cols)

  # Create Leaf based on Species + AccNum
  blast_out <- blast_out %>%
    mutate(sacc=paste0(sacc, ".1")) %>%
    select(AccNum=sacc, Species=sscinames) %>%
    separate(Species, into=c("Genus", "Spp"),
           sep=" ", remove=F,
           extra = "merge", fill = "left") %>%
    # 1 char from Genus, 3 char from species
    # kingdomPhylum_GenusSpecies
    mutate(Leaf=paste(paste0(str_sub(Genus, start=1, end=1),
                             str_sub(Spp, start=1, end=3)),
                      AccNum,
                      sep="_"))
  # blast_out$Leaf <- map(blast_out$Leaf, to_titlecase) # this returns a list class var

  #colnames(ipr_out) <- ipr_cols
  ipr_out$Strand <- rep("forward", nrow(ipr_out))

  ## GET AccNum GeneName, Species, Lineage mapping from the BLAST file!
  ipr_out <- ipr_out %>% left_join(blast_out, by="AccNum")

  ipr_out %>%
    group_by(Analysis) %>% summarize(occurrence=n()) %>%
    arrange(-occurrence)

  ipr_out <- ipr_out %>% arrange(Leaf, StartLoc, StopLoc)
  ipr_out_sub <- filter(ipr_out, grepl(pattern=analysis, x=Analysis))

  ## PLOTTING
  ## domains as separate arrows
  ggplot(ipr_out_sub,
         aes(xmin = StartLoc, xmax = StopLoc,
             y = Leaf, fill = SignDesc, label=SignAcc)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                    arrowhead_width = unit(1, "mm")) +
    geom_gene_label(align = "left") +
    # geom_blank(data = dummies) +
    facet_wrap(~ Analysis) + #, ncol = 1 + #scales = "free",
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position="bottom") +
    theme_minimal() + theme_genes()

  ## DA within genes
  ggplot(ipr_out_sub, fill="beige",
         #subset(ipr_out, molecule == "Genome4" & AccNum == "YP_177903.1"),
         aes(xmin = 1, xmax = SeqLen, y = Analysis)) +
    geom_gene_arrow() +
    # geom_gene_label(aes(label = Leaf)) +
    geom_subgene_arrow(
      data = ipr_out_sub, #molecule == "Genome4" & AccNum == "YP_177903.1"),
      aes(xsubmin = StartLoc, xsubmax = StopLoc, fill = SignDesc)
    ) +
    geom_subgene_label(
      data = ipr_out_sub, #molecule == "Genome4" & AccNum == "YP_177903.1"),
      aes(xsubmin = StartLoc, xsubmax = StopLoc, label = SignAcc),
      min.size = 0) +
    scale_fill_brewer(palette = "Set3") +
    facet_grid(.~Leaf, labeller=label_both, scales='free_x') +
    theme_minimal() + theme_genes()
}

##################################
## Description of gggenes colnames
##################################
## ggenes package:
## https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html
## Column Names FROM gggenes package for reference + mapping
## `example_genes`
## example: Genome5  genA 405113 407035 forward         1
#gggenes_cols <- c("molecule", "gene", "start", "end",
#                  "strand", "direction")
## `example_subgenes`
## example: Genome5  genA 405113 407035 forward  genA-1 405774 406538
#gggenes_sub_cols <- c("molecule", "gene", "start", "end",
#                      "strand", "subgene", "from", "to")

###############################
## DESCRIPTION of `ipr_cols` ##
###############################
## FOR INTERPROSCAN tsv files
## REF: https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats

# Protein Accession (e.g. P51587)
# Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
# Sequence Length (e.g. 3418)
# Analysis (e.g. Pfam / PRINTS / Gene3D)
# Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
# Signature Description (e.g. BRCA2 repeat profile)
# Start location
# Stop location
# Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
# Status - is the status of the match (T: true)
# Date - is the date of the run
# (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
# (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
# (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
# (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
