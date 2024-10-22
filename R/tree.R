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

##############################
## Approach 0-4
############
#' Approach 0: !! FastTree will only work if there are unique sequence names!!
#' createFA2Tree
#' @author Klangina, Outreachy Applicant
#' @author Janani Ravi, MolEcologist
#' @keywords phylogenetic tree, alignment, fasta, FastTree
#' @description
#' Generates a phylogenetic tree from an alignment file '.fa' using either FastTree or a detailed phylogenetic analysis approach.
#'
#' @param fa_file Character. Path to the alignment FASTA file (.fa) from which
#'  the phylogenetic tree will be generated. Default is "data/alns/pspa_snf7.fa".
#' @param out_file Character. Path to the output file where the generated tree (.tre) will 
#' be saved. Default is "data/alns/pspa_snf7.tre".
#' @param enable_fast Logical. If TRUE, uses FastTree approach; if FALSE, uses detailed phylogenetic analysis. Default is TRUE.
#' @param fasttree_path Character. Path to the FastTree executable. Default is "src/FastTree".
#' @param format Character. Format of the input file. Default is "fasta".
#' @param seq_type Character. Type of sequences in the input file. Default is "AA" (amino acid).
#' @param dist_matrix Character. Distance matrix to be used. Default is "similarity".
#' @param subset_size Numeric. Number of sequences to subset for analysis. Default is 10.
#' @param model_test Character. Model to use for model testing. Default is "all".
#' @param dna_model Character. Model to use for DNA distance calculation. Default is "JC69".
#' @param ml_model Character. Model to use for maximum likelihood estimation. Default is "JC".
#' @param rearrangement Character. Type of tree rearrangement for ML optimization. Default is "stochastic".
#' @param bootstrap_reps Numeric. Number of bootstrap replicates. Default is 100.
#' @param bootstrap_p Numeric. P-value threshold for bootstrap support. Default is 50.
#' @param plot_type Character. Type of plot for bootstrap tree. Default is "p".
#' @param bootstrap_file Character. Path to save the bootstrap tree file. Default is "bootstrap_example.tre".
#'
#' @importFrom ape write.tree
#' @importFrom phangorn bootstrap.pml dist.ml NJ modelTest phyDat plotBS pml pml.control pratchet optim.parsimony optim.pml read.phyDat upgma
#' @importFrom seqinr dist.alignment read.alignment
#' @importFrom stats logLik
#'
#' @return No return value. The function generates phylogenetic tree files 
#' (.tre) based on either FastTree or detailed phylogenetic analysis approaches. 
#' It also produces various plots and console outputs when using the detailed approach.
#' @export
#'
#' @details The function offers two main approaches for phylogenetic tree construction:
#' 1. FastTree approach: A quick method using the FastTree program.
#' 2. Detailed phylogenetic analysis: Includes model testing, distance matrix calculation,
#'    UPGMA and NJ tree construction, parsimony searches, and maximum likelihood estimation
#'    with bootstrapping.
#' The alignment file should be in FASTA format with sequences properly aligned.
#' When using FastTree, ensure that the FastTree executable is correctly installed and 
#' the path is properly specified.
#'
#' @note This function requires the ape, phangorn, and seqinr packages to be installed.
#' The detailed analysis approach performs computationally intensive operations and may 
#' take a significant amount of time to run, especially for large datasets or with high 
#' numbers of bootstrap replicates. FastTree approach is significantly faster but may 
#' provide less detailed results.
#'
#' @examples
#' \dontrun{
#' # Using FastTree approach
#' createFA2Tree("data/alns/my_alignment.fa", "data/trees/my_tree.tre", enable_fast = TRUE)
#' 
#' # Using detailed phylogenetic analysis approach
#' createFA2Tree("data/alns/protein_align.fa", "data/trees/protein_tree.tre",
#'               enable_fast = FALSE, subset_size = 20, bootstrap_reps = 500)
#' }
createFA2Tree <- function(
    fa_file = "data/alns/pspa_snf7.fa",
    out_file = "data/alns/pspa_snf7.tre",
    enable_fast = TRUE,
    fasttree_path = "src/FastTree",
    format = "fasta",
    seq_type = "AA",
    dist_matrix = "similarity",
    subset_size = 10,
    model_test = "all",
    dna_model = "JC69",
    ml_model = "JC",
    rearrangement = "stochastic",
    bootstrap_reps = 100,
    bootstrap_p = 50,
    plot_type = "p",
    bootstrap_file = "bootstrap_example.tre"
) {
    ## SAMPLE ARGS
    # fa_file="data/alns/pspa_snf7.fa"
    # out_file="data/alns/pspa_snf7.tre"
    # fasttree_path="src/FastTree"

    if (enable_fast) {
        # FastTree approach
        print(fa_file)
        system2(
            command = fasttree_path,
            args = paste(c(fa_file, ">", out_file),
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
    } else {
        # Detailed phylogenetic analysis approach
        ###########################
        ## Approach 1
        ###########################
        ## https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
        prot_aln <- read.alignment(file = fa_file, format = format)
        prot_aln$seq
        prot_dist <- dist.alignment(prot_aln, matrix = dist_matrix)
        prot_nj <- NJ(prot_dist)
        prot_nj_tree <- plot(prot_nj, main = "Neighbor Joining")
        write.tree(phy = prot_nj_tree, file = out_file)

        ###########################
        ## Approach 2
        ###########################
        ## Alignment file formats and conversion
        # read in sequence data, convert to phyDat
        prot_fa <- read.phyDat(fa_file, format = format, type = seq_type)
        prot_phyDat <- phyDat(prot_fa, type = seq_type, levels = NULL)
        prot_subset <- subset(prot_phyDat, 1:subset_size)
        prot_subset_phyDat <- phyDat(prot_subset, type = seq_type, levels = NULL)

        ## Model Testing & Distance Matrices
        ## Comparison of different nucleotide or amino acid substitution models
        mt <- modelTest(prot_subset, model = model_test)
        print(mt)

        # estimate a distance matrix using a Jules-Cantor Model
        dna_dist <- dist.ml(prot_subset, model = dna_model)

        ## Neighbor Joining, UPGMA, and Maximum Parsimony
        # quick and dirty UPGMA and NJ trees
        prot_UPGMA <- upgma(dna_dist)
        prot_NJ <- NJ(dna_dist)
        plot(prot_UPGMA, main = "UPGMA")
        plot(prot_NJ, main = "Neighbor Joining")

        # parsimony searches
        prot_optim <- optim.parsimony(prot_NJ, prot_phyDat)
        prot_pratchet <- pratchet(prot_subset) # returning error
        plot(prot_optim)
        plot(prot_pratchet)

        ## Maximum likelihood and Bootstrapping
        # ml estimation w/ distance matrix
        fit <- pml(prot_NJ, prot_subset)
        print(fit)
        fitJC <- optim.pml(fit, model = ml_model, rearrangement = rearrangement)
        logLik(fitJC)
        bs <- bootstrap.pml(fitJC,
            bs = bootstrap_reps, optNni = TRUE, multicore = TRUE,
            control = pml.control(trace = 0)
        )
        plotBS(midpoint(fitJC$tree), bs, p = bootstrap_p, type = plot_type)

        # subsetted alignment bs example
        prot_subset_dm <- dist.ml(prot_subset)
        prot_subset_NJ <- NJ(prot_subset_dm)
        fit2 <- pml(prot_subset_NJ, data = prot_subset)
        print(fit2)
        fitJC2 <- optim.pml(fit2, model = ml_model, rearrangement = rearrangement)
        logLik(fitJC2)
        bs_subset <- bootstrap.pml(fitJC2,
            bs = bootstrap_reps, optNni = TRUE, multicore = TRUE,
            control = pml.control(trace = 0)
        )
        bs2 <- plotBS(midpoint(fitJC2$tree), bs, p = bootstrap_p, type = plot_type)

        ## Exporting Trees
        write.tree(bs2, file = bootstrap_file)
    }
}
#' convertAlignment2Trees
#' @author Klangina, Outreachy Applicant
#' @description
#' Generate Trees for ALL fasta files in the specified directory. This function processes
#' multiple alignment files and creates phylogenetic trees for each, using either FastTree
#' or a detailed phylogenetic analysis approach.
#'
#' @param aln_path Character. Path to the directory containing all the alignment FASTA 
#' files (.fa) for which trees will be generated. Default is here("data/alns/").
#' @param enable_fast Logical. If TRUE, uses FastTree approach; if FALSE, uses detailed 
#' phylogenetic analysis. Default is TRUE.
#' @param fasttree_path Character. Path to the FastTree executable. Default is here("src/FastTree").
#' @param format Character. Format of the input files. Default is "fasta".
#' @param seq_type Character. Type of sequences in the input files. Default is "AA" (amino acid).
#' @param dist_matrix Character. Distance matrix to be used. Default is "similarity".
#' @param subset_size Numeric. Number of sequences to subset for analysis. Default is 10.
#' @param model_test Character. Model to use for model testing. Default is "all".
#' @param dna_model Character. Model to use for DNA distance calculation. Default is "JC69".
#' @param ml_model Character. Model to use for maximum likelihood estimation. Default is "JC".
#' @param rearrangement Character. Type of tree rearrangement for ML optimization. Default is "stochastic".
#' @param bootstrap_reps Numeric. Number of bootstrap replicates. Default is 100.
#' @param bootstrap_p Numeric. P-value threshold for bootstrap support. Default is 50.
#' @param plot_type Character. Type of plot for bootstrap tree. Default is "p".
#' @param bootstrap_file_suffix Character. Suffix for bootstrap tree files. Default is "_bootstrap.tre".
#'
#' @importFrom here here
#' @importFrom purrr pmap
#' @importFrom stringr str_replace_all
#'
#' @return No return value. The function generates tree files (.tre) and bootstrap files 
#' for each alignment file in the specified directory.
#' @export
#'
#' @details This function processes all .fa files in the specified directory. For each file, 
#' it generates a phylogenetic tree using either FastTree (if enable_fast is TRUE) or a 
#' detailed phylogenetic analysis approach (if enable_fast is FALSE). The detailed approach 
#' includes model testing, distance matrix calculation, UPGMA and NJ tree construction, 
#' parsimony searches, and maximum likelihood estimation with bootstrapping.
#'
#' @note This function requires the ape, phangorn, and seqinr packages to be installed.
#' The detailed analysis approach performs computationally intensive operations and may 
#' take a significant amount of time to run, especially for large datasets or with high 
#' numbers of bootstrap replicates. FastTree approach is significantly faster but may 
#' provide less detailed results.
#'
#' @examples
#' \dontrun{
#' # Using FastTree approach for all alignments in the default directory
#' convertAlignment2Trees()
#' 
#' # Using detailed phylogenetic analysis for alignments in a specific directory
#' convertAlignment2Trees(aln_path = "path/to/alignments/", 
#'                        enable_fast = FALSE, 
#'                        subset_size = 20, 
#'                        bootstrap_reps = 500)
#' }
convertAlignment2Trees <- function(
    aln_path = here("data/alns/"),
    enable_fast = TRUE,
    fasttree_path = here("src/FastTree"),
    format = "fasta",
    seq_type = "AA",
    dist_matrix = "similarity",
    subset_size = 10,
    model_test = "all",
    dna_model = "JC69",
    ml_model = "JC",
    rearrangement = "stochastic",
    bootstrap_reps = 100,
    bootstrap_p = 50,
    plot_type = "p",
    bootstrap_file_suffix = "_bootstrap.tre"
) {
    # finding all fasta alignment files
    fa_filenames <- list.files(path = aln_path, pattern = "*.fa")
    fa_paths <- paste0(aln_path, fa_filenames)
    variable <- str_replace_all(basename(fa_filenames),
        pattern = ".fa", replacement = ""
    )

    ## Using purrr::pmap to generate trees for each of the fa files!
    fa2tre_args <- list(
        fa_file = fa_paths,
        out_file = paste0(aln_path, variable, ".tre"),
        enable_fast = enable_fast,
        fasttree_path = fasttree_path,
        format = format,
        seq_type = seq_type,
        dist_matrix = dist_matrix,
        subset_size = subset_size,
        model_test = model_test,
        dna_model = dna_model,
        ml_model = ml_model,
        rearrangement = rearrangement,
        bootstrap_reps = bootstrap_reps,
        bootstrap_p = bootstrap_p,
        plot_type = plot_type,
        bootstrap_file = paste0(aln_path, variable, bootstrap_file_suffix)
    )
    
    pmap(
        .l = fa2tre_args, 
        .f = createFA2Tree
    )
}