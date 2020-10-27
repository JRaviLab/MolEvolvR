install.packages(c("sourcetools", "geometry", "crosstalk", "vdiffr", "rnaturalearth", "DT", "devtools", "diffobj", "freetypeharfbuzz", "gdtools", "shiny", "testthat", "usethis"), dependencies=T)
install.packages(c("manipulateWidget", "quadprog", "fastmatch", "maps", "gtools", "mnormt", "covr", "geiger", "rgl"), dependencies=T)
install.packages(c("phytools", "phangorn"), dependencies=T)

## Other packages
install.packages(c("conflicted","here","tidyverse","pdftools","latexpdf","UpSetR","wordcloud", "rlang"), dependencies=T)
install.packages(c("gggenes", "ggraph", "igraph", "visNetwork", "ggdendro", "gganimate", "ggnetwork", "networkD3", "geomnet", "ggrepel", "cowplot", "ggtree"), dependencies=T)

install.packages(c("shinyBS", "pins", "shinyjs", "sunburstR", "sodium", "rentrez"))

install.packages(c("viridis", "msa", "furrr"))

{if(Sys.info()["sysname"] == "Windows"){
  install.packages("XML", type =  "win.binary")
}
else{
  install.packages("XML")
}}
install.packages("rentrez")


remotes::install_github("paulc91/shinyauthr")
remotes::install_github("mhahsler/rMSA", force=T)
