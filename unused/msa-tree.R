## Multiple Sequence alignment & Phylogenetic trees
## Created: Apr 08, 2019
## Janani Ravi (@jananiravi)
## Tested with PSP data | PspN-DUF3046 fasta files

## msa package: https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
## ggtree: https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
## ape: http://www.phytools.org/eqg/Exercise_3.2/
## phylogram: https://cran.r-project.org/web/packages/phylogram/vignettes/phylogram-vignette.html
## tidytree: https://yulab-smu.github.io/treedata-book/chapter2.html
## ggtree + MSA + heatmap: https://yulab-smu.github.io/treedata-book/chapter7.html
## http://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/

## Installing and loading packages
## install CRAN Task View for phylogenetics install.packages('ctv')
## library('ctv') install.views('Phylogenetics')
## BiocManager::install("msa")
## update.views('Phylogenetics')
## BiocManager::install("ggtree", version = "3.8")

## install.packages("phylogram")
## install.packages("dendextend")
## install.packages("tidytree")


## Loading packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(msa))
# library(seqinr)
suppressPackageStartupMessages(library(Biostrings))
# library(phylogram)
# library(phytools) # depends on phangorn which depends on quadprog. so!!!
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidytree))

## Input files: Fasta format
# my_seqs_file <- read_tsv("data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt")
# colnames(my_seqs_file) <- c("name", "sequence")
# my_alignment <- as(my_seqs_file, "BStringSet")
fasta_filepath <- "data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt"
my_seqs <- readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
my_seqs_msa <- msa(my_seqs)

## Print alignment to screen
print(my_seqs, show="complete")
## Print colorufl alignment to file
msaPrettyPrint(my_seqs_msa, output="pdf",
               file=paste0(fasta_filepath, ".pdf"),
               y=c(260,325),
               showNames="left", showLogo="top",
               logoColors="rasmol", # “chemical”, “rasmol”, “hydropathy”, “structure”, “standard area”, “accessible area”
               shadingMode="functional", # or "similar"
               shadingModeArg="structure",
               shadingColors="blues",
               consensusColors="ColdHot",
               askForOverwrite=FALSE, verbose=FALSE,
               furtherCode=c("\\defconsensus{.}{lower}{upper}",
                             "\\showruler{1}{top}"))

## using package seqinr
my_seqs_msa_aln <- msaConvert(my_seqs_msa, type="seqinr::alignment")

## From seqinr
# d <- dist.alignment(my_seqs_msa_aln, "identity")
# as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

## Phylogenetic tree
## using package ape
## build neighbor-joining tree
seqTree <- nj(d)
plot(seqTree, main="Phylogenetic Tree of MSA")

## drawing trees using ggtree
# ggtree(seqTree)
groupInfo <- split(seqTree$tip.label,
                   gsub("_\\w+", "", seqTree$tip.label))
seqTree <- groupOTU(seqTree, groupInfo)

ggtree(seqTree, aes(color=group),
       layout='circular') +
  geom_tiplab(size=1, aes(angle=angle))


ggtree(seqTree) + #xlim(NA, 6) +
  geom_cladelabel(node=1,
                  label="Mtb", align=T, color='red') +
  geom_cladelabel(node=32,
                  label="Osp", align=T, color='blue')


## Tree + MSA
# beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
# beast_tree <- read.beast(beast_file)
# fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
# msaplot(ggtree(beast_tree), fasta)

offs <- 0
msaplot(ggtree(seqTree), fasta=fasta_filepath, offset=0.75) +
  geom_tiplab(size=2, align=TRUE, linesize=.5) #+
  # geom_cladelabel(node=2,
  #                 label="Mtb", align=T, color='red', offset=offs) +
  # geom_cladelabel(node=12,
  #                 label="Msmeg", align=T, color='blue', offset=offs) +
  # geom_cladelabel(node=23,
  #                 label="Cptb", align=T, color='darkgreen', offset=offs)


##########################
#### Doesn't work yet ####
##########################
# ## root with XXX as outgroup
# phy1 <- root(seqTree, "A*.sp.FB24|Arthrobacter sp FB24")
# phy2 <- root(seqTree2, "XXX")
# ## convert phylo objects to dendrograms
# dnd1 <- as.cladogram(phy1)
# dnd2 <- as.dendrogram(phy2)
# ## rearrange in ladderized fashion
# dnd1 <- ladder(dnd1)
# dnd2 <- ladder(dnd2)
# ## plot the tanglegram
# dndlist <- dendextend::dendlist(dnd1, dnd2)
# dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)

#### Doesn't work either!!
# ## trying to root with H37Rv
# seqTree.dend <- as.cladogram(seqTree)
# ## isolate root node (species H37Rv)
# ancestor <- prune(seqTree.dend, pattern = "Arthrobacter", keep = TRUE)
# ## alternative option using subset operator
# # ancestor <- x[[2]][[2]]
# ## create subtree without species C
# subtree <- prune(seqTree.dend, pattern = "Arthrobacter", keep=FALSE)
# ## graft subtree onto root
# x <- list(ancestor, subtree)
# ## set attributes as above
# x <- as.cladogram(remidpoint(x))
# ## plot dendrogram
# plot(x, yaxt = "n")
