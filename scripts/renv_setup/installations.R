
# CRAN installations
install.packages(c('devtools', "usethis", "testthat", "reprex",
                   "docstring", "roxygen2"),
                 dependencies=T)

install.packages(c("conflicted", "tidyverse", "here",
                   "rmarkdown", "knitr", "xtable", "DT", "gt",
                   "blogdown", "bookdown", "pkgdown",
                   "shiny", "shinydashboard", "flexdashboard"),
                 dependencies=T)

install.packages(c("Rcpp", "data.table", "rlang", "glue",
                   "future", "furrr",
                   "remotes", "rgl",
                   "htmlwidgets", "gh", "httr", "XML",
                   "drat", "rhub", "covr",
                   "RPushbullet", "jsonlite"),
                 dependencies=T)

install.packages(c("gridExtra", "cowplot", "gganimate", "plotly", "tidytext",
                   "ggvis", "ggthemes", "ggsci", "gggenes",
                   "UpSetR", "wordcloud2", "wordcloud", "sunburstR",
                   "igraph", "ggraph", "visNetwork", #"network3D", "d3Network",
                   "heatmap3", "heatmaply", "viridis", "d3r",
                   "rentrez", "reutils", "biomartr",
                   "ape", "phylogram", "phangorn", "seqinr",
                   "phylotools", "phytools", "tidytree",
                   "pdftools", "latexpdf", "tinytex"),
                 dependencies=T)

# Bioconductor installations
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies=T)
BiocManager::install("msa")
BiocManager::install("Biostrings")

# GitHub and Devtools installations
devtools::install_github("rstudio/d3heatmap")
remotes::install_github("mhahsler/rMSA")
devtools::install_github("GuangchuangYu/ggtree")

# Other
remotes::install_github("ropenscilabs/travis")

# trouble with OpenMP
# BiocManager::install("ggtree") # ERROR w/ fopenmp