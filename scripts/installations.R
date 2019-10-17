
# CRAN installations
install.packages(c("devtools", "BiocManager",
                   "conflicted","docstring", "roxygen2",
                   "tidyverse", "here", "testthat", "reprex", "data.table",
                   "rmarkdown", "knitr", "blogdown", "bookdown",
                   "Rcpp", "rgl",
                   "gridExtra", "cowplot", "gganimate", "plotly", "ggthemes",
                   "UpSetR", "wordcloud2", "igraph", "tidytext",
                   "ape", "phylogram", "msa", "seqinr", "phylotools",
                   "tidytree", "ggtree"), dependencies=T)

# Bioconductor installations
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# GitHub and Devtools installations
