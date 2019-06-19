library(ape)
library(ggtree)
library(tidytree)
library(seqinr)

#take in a file, generate trees from that
seq_tree <- function(fasta_filepath){
  my_seqs <- readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
  my_seqs_msa <- msa(my_seqs)
  my_seqs_msa_aln <- msaConvert(my_seqs_msa, type="seqinr::alignment")

  #below was commented out, does it need to change as one of the parameters? the bottom keeps
  d <- dist.alignment(my_seqs_msa_aln, "identity")
  #as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

  ## Phylogenetic tree
  ## using package ape
  ## build neighbor-joining tree
  seqTree <- nj(d)
  #plot(seqTree, main="Phylogenetic Tree of MSA")

  groupInfo <- split(seqTree$tip.label,
                     gsub("_\\w+", "", seqTree$tip.label))
  seqTree <- groupOTU(seqTree, groupInfo)

  # ggtree(seqTree, aes(color=group),
  #        layout='circular') +
  #   geom_tiplab(size=1, aes(angle=angle))

  offs <- 0
  msaplot(ggtree(seqTree), fasta=fasta_filepath, offset=0.75) +
    geom_tiplab(size=2, align=TRUE, linesize=.5)

}

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