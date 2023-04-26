suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))

## REFS
## https://sashamaps.net/docs/resources/20-colors/
## https://mokole.com/palette.html # (distinct colors)
## https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html
## http://www.sthda.com/english/wiki/colors-in-r
## https://github.com/karthik/wesanderson

# discrete_palettes <- c( # arc>bac>euk>vir
#   RColorBrewer::brewer.pal(8, "Set2"),                         # 1-8
#   'darksalmon', 'rosybrown', 'navy', 'firebrick3',             # 9-12
#   'yellowgreen', 'deepskyblue3', 'orange', 'gold3',            # 13-16
#   'royalblue1', 'cyan4', 'darkmagenta', 'khaki',               # 18-20
#   'salmon', 'darkolivegreen3', 'powderblue', 'mediumpurple1',  # 21-24
#   RColorBrewer::brewer.pal(7, "Dark2")                         # 25-31
# )
#
#
# lineages_short <- c("Vir", "Euk", "Bacteria", "B.thermot",
#              "B.thermob", "B.teneri", "B.spiro", "B.PVC",
#              "B.proteo", "B.nitro", "B.igna", "B.gemma",
#              "B.fuso", "B.firmi", "B.FCB", "B.elusi",
#              "B.dictyo", "B.deino.therm", "B.cyano" , "B.chloroflexi",
#              "B.chlorobi", "B.caldit", "B.caldis", "B.arma",
#              "B.aqui", "B.actino", "B.acido", "Archaea",
#              "A.TACK", "A.eury", "A.asgard")
#
# lineages <- c("Viruses", "Eukaryota", "Bacteria",
#               "Bac>thermotogae", "Bac>thermobaculum", "Bac>tenericutes",
#               "Bac>spirochaetes", "Bac>PVC_group", "Bac>proteobacteria",
#               "Bac>nitrospirae", "Bac>ignavibacteriae", "Bac>gemmatimonadetes",
#               "Bac>fusobacteria", "Bac>firmicutes", "Bac>FCB_group",
#               "Bac>elusimicrobia", "Bac>dictyoglomi", "Bac>deinococcus-thermus",
#               "Bac>cyanobacteria", "Bac>chloroflexi", "Bac>chlorobi",
#               "Bac>calditrichaeota", "Bac>caldiserica", "Bac>armatimonadetes",
#               "Bac>aquificae", "Bac>actinobacteria", "Bac>acidobacteria",
#               "Archaea", "Arc>TACK_group", "Arc>euryarchaeota",
#               "Arc>asgard_group")
#
# lin_color <- tibble(lin=lineages, col=discrete_palettes)

lin_color <- read_tsv("data/acc_files/lin_color.tsv")
par(mar = rep(0, 4))
pie(rep(1, length(lin_color$col)), init.angle=50,
    col=lin_color$col,
    labels=rev(lin_color$lin))
