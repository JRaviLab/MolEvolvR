## Generating Phylogenetic Trees from Alignment Fasta files
## Includes the following functions:
## convertAlignment2Trees, convertFA2Tree, createFA2Tree
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
#' convertFA2Tree
#'
#' @param fa_path Path to the input FASTA alignment file (.fa). Default is the 
#' path to "data/alns/pspa_snf7.fa".
#' @param tre_path Path to the output file where the generated tree (.tre) will 
#' be saved. Default is the path to "data/alns/pspa_snf7.tre".
#' @param fasttree_path Path to the FastTree executable, which is used to 
#' generate the phylogenetic tree. Default is "src/FastTree".
#' 
#' @importFrom rlang abort
#'
#' @return No return value. The function generates a tree file (.tre) from the 
#' input FASTA file.
#' @export
#'
#' @examples
#' \dontrun{
#' convert_fa2tre(here("data/alns/pspa_snf7.fa"), 
#'                 here("data/alns/pspa_snf7.tre"), 
#'                 here("src/FastTree")
#' }
convertFA2Tree <- function(fa_path = here("data/alns/pspa_snf7.fa"),
    tre_path = here("data/alns/pspa_snf7.tre"),
    fasttree_path = here("src/FastTree")) {
    # fa_path=here("data/alns/pspa_snf7.fa")
    # tre_path=here("data/alns/pspa_snf7.tre")
    # fasttree_path=here("src/FastTree")
    
    # Check if the FASTA file exists
    if (!file.exists(fa_path)) {
        abort(paste("Error: The FASTA file does not exist at:", fa_path))
    }
    
    # Check if the FastTree executable exists
    if (!file.exists(fasttree_path)) {
        abort(paste("Error: The FastTree executable does not exist at:", 
                   fasttree_path))
    }
    
    # Check if the output directory exists
    tre_dir <- dirname(tre_path)
    if (!dir.exists(tre_dir)) {
        abort(paste("Error: The output directory does not exist:", tre_dir))
    }
    
    # Check if the output file already exists
    if (file.exists(tre_path)) {
        cat("Warning: The output file already exists and will be overwritten:", 
            tre_path, "\n")
    }
    
    message(fa_path)
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
#' convertAlignment2Trees
#'
#' @description
#' Generate Trees for ALL fasta files in "data/alns"
#'
#' @param aln_path Path to the directory containing all the alignment FASTA 
#' files (.fa) for which trees will be generated. Default is "data/alns/".
#' 
#'
#' @importFrom here here
#' @importFrom purrr pmap
#' @importFrom stringr str_replace_all
#'
#' @return No return value. The function generates tree files (.tre) for each 
#' alignment file in the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_trees(here("data/alns/"))
#' }
convertAlignment2Trees <- function(aln_path = here("data/alns/")) {
    
    # Check if the alignment directory exists
    if (!dir.exists(aln_path)) {
        abort(paste("Error: The alignment directory does not exist:", aln_path))
    }
    # finding all fasta alignment files
    fa_filenames <- list.files(path = aln_path, pattern = "*.fa")
    # Check if any FASTA files were found
    if (length(fa_filenames) == 0) {
        abort("Error: No FASTA files found in the specified directory.")
    }
    
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
        .l = fa2tre_args, .f = convertFA2Tree,
        fasttree_path = here("src/FastTree")
    )
}

##############################
## REFS: 1-4
############
#' createFA2Tree
#'
#' @author Janani Ravi, MolEcologist
#' @keywords phylogenetic tree, alignment, fasta
#' @description
#' Generating phylogenetic tree from alignment file '.fa'
#'
#' @param fa_file Character. Path to the alignment FASTA file (.fa) from which
#'  the phylogenetic tree will be generated. Default is 'pspa_snf7.fa'.
#' @param out_file Path to the output file where the generated tree (.tre) will 
#' be saved. Default is "data/alns/pspa_snf7.tre".
#'
#' @importFrom ape write.tree
#' @importFrom phangorn bootstrap.pml dist.ml NJ modelTest phyDat plotBS pml pml.control pratchet optim.parsimony optim.pml read.phyDat upgma
#' @importFrom seqinr dist.alignment read.alignment
#' @importFrom stats logLik
#'
#' @return No return value. The function generates a phylogenetic tree file 
#' (.tre) based on different approaches like Neighbor Joining, UPGMA, and 
#' Maximum Likelihood.
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
createFA2Tree <- function(fa_file = "data/alns/pspa_snf7.fa",
    out_file = "data/alns/pspa_snf7.tre") {
    ## SAMPLE ARGS
    # fa_file="data/alns/pspa_snf7.fa"
    # out_file="data/alns/pspa_snf7.tre"
    
    # Check if the FASTA file exists
    if (!file.exists(fa_file)) {
        abort(paste("Error: The FASTA file does not exist at:", fa_file))
    }
    
    # Check if the output directory exists
    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) {
        abort(paste("Error: The output directory does not exist:", out_dir))
    }
    
    # Check if the output file already exists
    if (file.exists(out_file)) {
        cat("Warning: The output file already exists and will be overwritten:", 
            out_file, "\n")
    }

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
    message(mt)

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
    message(fit)
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
    message(fit2)
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
