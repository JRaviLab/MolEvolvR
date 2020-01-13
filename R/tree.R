## Generate Phylogenetic Tree
## Includes the following functions:
## generate_fa2tre
## Modified: Dec 27, 2019
## Molecular Ecologist (@molecologist), Janani Ravi (@jananiravi)

#################
## Pkgs needed ##
#################
library(ape)
library(phangorn)
library(seqinr)
library(here); library(tidyverse)
library(data.table)
conflicted::conflict_prefer("filter", "dplyr")

###############
## References #
###############
## 1a. Quick & Dirty Tree building in R: https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/ &
## 1b: Script: https://github.com/elinck/molecologist/blob/master/molecoltutorial_treesinR.R
## 2. ape: https://cran.r-project.org/web/packages/ape/ape.pdf
## 3a. phangorn: https://cran.r-project.org/web/packages/phangorn/phangorn.pdf
## 3b. phangorn vignette: https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
## 4. MSA & Phylogenetic Trees (NJ, rooted, unrooted, retrieve sequences, ... ) https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html

###############
## Functions ##
###############
generate_fa2tre <- function(fa_file = "data/alns/pspa_snf7.fa",
                            out_file = "data/alns/pspa_snf7.tre") {
  #' Generating phylogenetic tree from alignment file '.fa'
  #' @author Janani Ravi, MolEcologist
  #' @keywords phylogenetic tree, alignment, fasta
  #' @description ...
  #' @param fa_file Character. Path to file.
  #'  Default is 'pspa_snf7.fa'
  #' @param xyz Boolean. If TRUE, a ...
  #' @examples generate_aln2tree('pspa_snf7.fa')
  #' @details The alignment file would need two columns: 1. accession +
  #' number and 2. alignment. The protein homolog accession to lineage mapping +
  #' file should have
  #'
  #' @note Please refer to the source code if you have alternate +
  #' file formats and/or column names.

  ## SAMPLE ARGS
  # fa_file = "data/alns/pspa_snf7.fa"
  # out_file = "data/alns/pspa_snf7.tre"

  ###########################
  ## Approach 1
  ###########################
  ## https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
  prot_aln <- read.alignment(fa_file, format="fasta")
  prot_aln$seq
  prot_dist <- dist.alignment(prot_aln, matrix="similarity")
  prot_nj <- NJ(prot_dist)
  prot_nj_tree <- plot(prot_nj, main = "Neighbor Joining")
  write.tree(prot_nj_tree, "pspa_snf7_nj.tree")

  ###########################
  ## Approach 2
  ###########################
  ## Alignment file formats and conversion
  # read in sequence data, convert to phyDat
  prot_fa <- read.aa(fa_file, format="fasta")
  prot_phyDat <- phyDat(prot_fa, type = "AA", levels = NULL)
  prot10 <- subset(prot_phyDat, 1:10)
  prot10_phyDat <- phyDat(prot10, type = "AA", levels = NULL)

  ## Model Testing & Distance Matrices
  ## Comparison of different nucleotide or amino acid substitution models
  mt <- modelTest(prot10, model="all")
  print(mt)

  #estimate a distance matrix using a Jules-Cantor Model
  dna_dist <- dist.ml(prot10, model="JC69")

  ## Neighbor Joining, UPGMA, and Maximum Parsimony
  #quick and dirty UPGMA and NJ trees
  prot_UPGMA <- upgma(dna_dist)
  prot_NJ  <- NJ(dna_dist)
  plot(prot_UPGMA, main="UPGMA")
  plot(prot_NJ, main = "Neighbor Joining")

  #parsimony searches
  prot_optim <- optim.parsimony(prot_NJ, prot)
  prot_pratchet <- pratchet(prot10) #returning error
  plot(prot_optim)
  plot(prot_pratchet)

  ## Maximum likelihood and Bootstrapping
  #ml estimation w/ distance matrix
  fit <- pml(prot_NJ, prot10)
  print(fit)
  fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
  logLik(fitJC)
  bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
  plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

  #subsetted alignment bs example
  prot10_dm <- dist.ml(prot10)
  prot10_NJ  <- NJ(prot10_dm)
  fit2 <- pml(prot10_NJ, data=prot10)
  print(fit2)
  fitJC2 <- optim.pml(fit2, model = "JC", rearrangement = "stochastic")
  logLik(fitJC2)
  bs_subset <- bootstrap.pml(fitJC2, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
  bs2 <- plotBS(midpoint(fitJC2$tree), bs, p = 50, type="p")

  ## Exporting Trees
  write.tree(bs2, file="bootstrap_example.tre")
}