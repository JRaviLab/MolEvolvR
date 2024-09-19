## Generating Phylogenetic Trees from Alignment Fasta files
## Includes the following functions:
## generate_trees, convert_fa2tre, generate_fa2tre
## Modified: Jan, 2020
## Janani Ravi (@jananiravi), Molecular Ecologist (@molecologist)

#################
## Pkgs needed ##
#################
# Only these here+tidyerse are needed for FastTree
# suppressPackageStartupMessages(library(here))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(ape))
# suppressPackageStartupMessages(library(phangorn))
# suppressPackageStartupMessages(library(seqinr))
# conflicted::conflict_prefer("filter", "dplyr")

###############
## References #
###############
## 0. Using FastTree 2.0: http://www.microbesonline.org/fasttree/
## 1a. Quick & Dirty Tree building in R: https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/ &
## 1b: Script: https://github.com/elinck/molecologist/blob/master/molecoltutorial_treesinR.R
## 2. ape: https://cran.r-project.org/web/packages/ape/ape.pdf
## 3a. phangorn: https://cran.r-project.org/web/packages/phangorn/phangorn.pdf
## 3b. phangorn vignette: https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
## 4. MSA & Phylogenetic Trees (NJ, rooted, unrooted, retrieve sequences, ... ) https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html

###############
## Functions ##
###############

###########################
## Approach 0 | FastTree2.0
###########################
## !! FastTree will only work if there are unique sequence names!!
#' convert_fa2tre
#'
#' @param fa_path
#' @param tre_path
#' @param fasttree_path
#'
#' @return
#' @export
#'
#' @examples
convert_fa2tre <- function(fa_path = here("data/alns/pspa_snf7.fa"),
    tre_path = here("data/alns/pspa_snf7.tre"),
    fasttree_path = here("src/FastTree")) {
    # fa_path=here("data/alns/pspa_snf7.fa")
    # tre_path=here("data/alns/pspa_snf7.tre")
    # fasttree_path=here("src/FastTree")
    print(fa_path)
    system2(
        command = fasttree_path,
        args = paste(c(fa_path, ">", tre_path),
            sep = "", collapse = " "
        )
    )
    ################################
    ## Compiling FastTree.c on a Mac
    ################################
    ## Other platforms? Check http://www.microbesonline.org/fasttree/
    # system(paste("gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o",
    #              here("src/FastTree"),
    #              here("src/FastTree.c"), "-lm", collapse=" "))
}
## Generate Trees for ALL fasta files in "data/alns"
#' generate_trees
#'
#' @description
#' Generate Trees for ALL fasta files in "data/alns"
#'
#' @param aln_path
#'
#' @importFrom here here
#' @importFrom purrr pmap
#' @importFrom stringr str_replace_all
#'
#' @return
#' @export
#'
#' @examples
generate_trees <- function(aln_path = here("data/alns/")) {
    # finding all fasta alignment files
    fa_filenames <- list.files(path = aln_path, pattern = "*.fa")
    fa_paths <- paste0(aln_path, fa_filenames)
    variable <- str_replace_all(basename(fa_filenames),
        pattern = ".fa", replacement = ""
    )

    ## Using purrr::pmap to generate trees for each of the fa files!
    fa2tre_args <- list(
        fa_path = fa_paths,
        tre_path = paste0(aln_path, variable, ".tre")
    )
    pmap(
        .l = fa2tre_args, .f = convert_fa2tre,
        fasttree_path = here("src/FastTree")
    )
}

##############################
## REFS: 1-4
############
#' generate_fa2tre
#'
#' @author Janani Ravi, MolEcologist
#' @keywords phylogenetic tree, alignment, fasta
#' @description
#' Generating phylogenetic tree from alignment file '.fa'
#'
#' @param fa_file Character. Path to file.
#'  Default is 'pspa_snf7.fa'
#' @param out_file
#'
#' @importFrom ape write.tree
#' @importFrom phangorn bootstrap.pml dist.ml NJ modelTest phyDat plotBS pml pml.control pratchet optim.parsimony optim.pml read.phyDat upgma
#' @importFrom seqinr dist.alignment read.alignment
#' @importFrom stats logLik
#'
#' @return
#' @export
#'
#' @details The alignment file would need two columns: 1. accession +
#' number and 2. alignment. The protein homolog accession to lineage mapping +
#' file should have
#'
#' @note Please refer to the source code if you have alternate +
#' file formats and/or column names.
#'
#' @examples
#' \dontrun{
#' generate_aln2tree("pspa_snf7.fa")
#' }
generate_fa2tre <- function(fa_file = "data/alns/pspa_snf7.fa",
    out_file = "data/alns/pspa_snf7.tre") {
    ## SAMPLE ARGS
    # fa_file="data/alns/pspa_snf7.fa"
    # out_file="data/alns/pspa_snf7.tre"

    ###########################
    ## Approach 1
    ###########################
    ## https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
    prot_aln <- read.alignment(file = fa_file, format = "fasta")
    prot_aln$seq
    prot_dist <- dist.alignment(prot_aln, matrix = "similarity")
    prot_nj <- NJ(prot_dist)
    prot_nj_tree <- plot(prot_nj, main = "Neighbor Joining")
    write.tree(phy = prot_nj_tree, file = tre_path)

    ###########################
    ## Approach 2
    ###########################
    ## Alignment file formats and conversion
    # read in sequence data, convert to phyDat
    prot_fa <- read.phyDat(fa_file, format = "fasta", type = "AA")
    prot_phyDat <- phyDat(prot_fa, type = "AA", levels = NULL)
    prot10 <- subset(prot_phyDat, 1:10)
    prot10_phyDat <- phyDat(prot10, type = "AA", levels = NULL)

    ## Model Testing & Distance Matrices
    ## Comparison of different nucleotide or amino acid substitution models
    mt <- modelTest(prot10, model = "all")
    print(mt)

    # estimate a distance matrix using a Jules-Cantor Model
    dna_dist <- dist.ml(prot10, model = "JC69")

    ## Neighbor Joining, UPGMA, and Maximum Parsimony
    # quick and dirty UPGMA and NJ trees
    prot_UPGMA <- upgma(dna_dist)
    prot_NJ <- NJ(dna_dist)
    plot(prot_UPGMA, main = "UPGMA")
    plot(prot_NJ, main = "Neighbor Joining")

    # parsimony searches
    prot_optim <- optim.parsimony(prot_NJ, prot)
    prot_pratchet <- pratchet(prot10) # returning error
    plot(prot_optim)
    plot(prot_pratchet)

    ## Maximum likelihood and Bootstrapping
    # ml estimation w/ distance matrix
    fit <- pml(prot_NJ, prot10)
    print(fit)
    fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
    logLik(fitJC)
    bs <- bootstrap.pml(fitJC,
        bs = 100, optNni = TRUE, multicore = TRUE,
        control = pml.control(trace = 0)
    )
    plotBS(midpoint(fitJC$tree), bs, p = 50, type = "p")

    # subsetted alignment bs example
    prot10_dm <- dist.ml(prot10)
    prot10_NJ <- NJ(prot10_dm)
    fit2 <- pml(prot10_NJ, data = prot10)
    print(fit2)
    fitJC2 <- optim.pml(fit2, model = "JC", rearrangement = "stochastic")
    logLik(fitJC2)
    bs_subset <- bootstrap.pml(fitJC2,
        bs = 100, optNni = TRUE, multicore = TRUE,
        control = pml.control(trace = 0)
    )
    bs2 <- plotBS(midpoint(fitJC2$tree), bs, p = 50, type = "p")

    ## Exporting Trees
    write.tree(bs2, file = "bootstrap_example.tre")
}
