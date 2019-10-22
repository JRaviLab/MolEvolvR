
# CRAN installations
install.packages(c("devtools", "conflicted","docstring", "roxygen2",
                   "tidyverse", "here", "testthat", "reprex",
                   "rmarkdown", "knitr", "blogdown", "bookdown",
                   "Rcpp", "rgl",
                   "gridExtra", "cowplot", "gganimate", "plotly",
                   "ggthemes", "ggsci",
                   "UpSetR", "wordcloud2", "igraph", "tidytext",
                   "ape", "phylogram", "seqinr", "phylotools", "phytools",
                   "tidytree"), dependencies=T)

install.packages("BiocManager", dependencies=T)
# trouble with OpenMP
# install.packages("data.table", dependencies=T)

# Bioconductor installations
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("msa")
# BiocManager::install("ggtree") # ERROR w/ fopenmp

# GitHub and Devtools installations
devtools::install_github("GuangchuangYu/ggtree")
