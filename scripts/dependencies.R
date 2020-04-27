install.packages(c("sourcetools", "geometry", "crosstalk", "vdiffr", "rnaturalearth", "DT", "devtools", "diffobj", "freetypeharfbuzz", "gdtools", "shiny", "testthat", "usethis"), dependencies=T)
install.packages(c("manipulateWidget", "quadprog", "fastmatch", "maps", "gtools", "mnormt", "covr", "geiger", "rgl"), dependencies=T)
install.packages(c("phytools", "phangorn"), dependencies=T)

## Other packages
install.packages(c("conflicted","here","tidyverse","pdftools","latexpdf","UpSetR","wordcloud"), dependencies=T)
install.packages(c("gggenes", "ggraph", "igraph", "ggdendro", "gganimate", "ggnetwork", "networkD3", "geomnet", "ggrepel", "cowplot", "ggtree"), dependencies=T)

remotes::install_github("mhahsler/rMSA", force=T)
