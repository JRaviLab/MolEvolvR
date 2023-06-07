## SCRIPT to Plot Phylogenetic Tree and MSA w/ protein FASTA files
## Created: May 13, 2019 | Janani Ravi
## Last Edited: May 14

###################
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidytree))
##################

sessionInfo()

## Input files: Fasta format
# my_seqs_file <- read_tsv("data/alignments/pspn-duf3046-aln/pspn.31seq.aln.txt")
# colnames(my_seqs_file) <- c("name", "sequence")
# my_alignment <- as(my_seqs_file, "BStringSet")
fasta_filepath <- "data/map3773c/all_seq_seqinr_format.fa"
my_seqs <- Biostrings::readAAStringSet(fasta_filepath) #, format="fasta", seek.first.rec=T)
my_seqs_msa <- msa::msa(my_seqs)

## Print alignment to screen
print(my_seqs, show="complete")
## Print colorufl alignment to file
msa::msaPrettyPrint(my_seqs_msa, output="pdf",
               file=paste0(fasta_filepath, ".pdf"),
               y=c(1,90),
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
my_seqs_msa_aln <- msa::msaConvert(my_seqs_msa, type="seqinr::alignment")

## From seqinr
d <- seqinr::dist.alignment(my_seqs_msa_aln, "identity")
# as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]

## Phylogenetic tree
## using package ape
## build neighbor-joining tree
seqTree <- ape::nj(d)
plot(seqTree, main="Phylogenetic Tree of MSA")

## drawing trees using ggtree
# ggtree(seqTree)
groupInfo <- split(seqTree$tip.label,
                   gsub("_\\w+", "", seqTree$tip.label))
seqTree <- tidytree::groupOTU(seqTree, groupInfo)


## Tree + MSA

offs <- 0
msaplot(ggtree(seqTree), fasta=fasta_filepath, offset=0.75) +
  geom_tiplab(size=2, align=F, linesize=.5)
