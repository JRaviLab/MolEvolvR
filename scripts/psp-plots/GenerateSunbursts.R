library(tidyverse)
library(sunburstR)
library(RColorBrewer)
library(d3r)
source("R/cleanup.R")

# read in all and clean
all <- read_tsv("data/rawdata_tsv/all_semiclean.txt")

colnames(all)[which(colnames(all) == "DomArch")] = "DomArch.ind"
colnames(all)[which(colnames(all) == "DomArch.orig")] = "DomArch.ind.orig"

colnames(all)[which(colnames(all) == "ClustName")] = "DomArch"

all <- all %>%
  cleanup_clust(domains_keep = domains_keep, domains_rename = domains_rename,repeat2s = FALSE)
colnames(all)[which(colnames(all) == "ClustName")] = "DomArch.repeats"

colnames(all)[which(colnames(all) == "ClustName.orig")] = "DomArch.orig"

colnames(all)[which(colnames(all) == "GenContext")] = "Temp"

all <- cleanup_gencontext(all,domains_rename = domains_rename, repeat2s = FALSE)
colnames(all)[which(colnames(all) == "GenContext")] = "GenContext.repeats"
colnames(all)[which(colnames(all) =="Temp")] = "GenContext"


all <- all %>% select(AccNum, DomArch, DomArch.repeats, GenContext, GenContext.repeats,
                      Lineage, Species, GeneName, Length, GCA_ID) %>% distinct()


# all requires columns Lineage and Domarch/Clustname




prot_all <- all

prot_all$Lineage = unlist(map(prot_all$Lineage, function(x){
  gt_pos = gregexpr(pattern = ">", x)[[1]][2]
  # If there is not a second '>' in the lineage
  if(is.na(gt_pos))
  {
    x
  }
  else{
    substr(x,1, (gt_pos-1) )
  }
} ))


df <- prot_all

rownames(df) <- NULL

legend_items <- unique(unlist(strsplit(df$Lineage, ">")))

## Assign Colors:
#cols <- sample(colorRampPalette(brewer.pal(9, 'Set1'))(length(legend_items)))
cols <- sample(colorRampPalette(brewer.pal(12, 'Paired'))(length(legend_items)))
# cols <- sample(colorRampPalette(brewer.pal(12, 'Pastel2'))(length(legend_items)))

queries = c("PspA", "Snf7", "PspB","PspC", "PspM",  "PspN", "LiaI-LiaF-TM", "Toast-rack")

# cols %>% paste0(collapse = "','")

## Current Colors
cols = c(
  'F48829','#ff4545','#7F9D55','#ff69e3','#D49B84',
  '#8963AD','#F6F399','#A3D77F','#AA9C6C','#6CA9CF',
  '#E62B2C','#1e00b3','#D29F57','#b3ff00','#754AA0',
  '#3385BB','#E2C26F','#A6CEE3','#F06161','#F58F57',
  '#EA4933','#977899','#D6CA99','#B15928','#F57C7C',
  '#5097C5','#37A22F','#34ebd2','#F3E587','#52AF43',
  '#A3D58E','#4693A8','#FDB660','#277DB1','#C6ADD3',
  '#D3A9B0','#65A99F','#84BF96','#FDA848','#FA9796',
  '#DE9E83','#784F99','#F06C45','#EB4647','#C17C3F',
  '#e28cff','#FE8D19','#9D7BBA','#B7A199','#559E3E',
  '#0dff00','#6DBD57','#FBB268')

# library(hash)
# find_dups <- function(c_vec)
# {
#   h <- hash()
#   for(i in c_vec)
#   {
#     if(has.key(i, h))
#     {
#       h[[i]] = h[[i]] +1
#     }
#     else
#     {
#       h[[i]] = 1
#     }
#   }
#   for(k in keys(h))
#   {
#     if( (h[[k]]) == 1)
#     {
#       del(k, h)
#     }
#   }
#   return(h)
# }


# df$Lineage <- str_replace_all(df$Lineage, "-", "_")
#
# df$relationship <- str_replace_all(df$Lineage, ">", "-")

# plot each lineage
plt_fun <- function(query){

  levels_vec = c("Level1", "Level2")

  # filter for dt
  plt_df <- df %>%
    filter(grepl(query, DomArch, ignore.case = T))

  # Take lineage column and break into the first to levels
  prot <- plt_df %>% select(Lineage)
  protLevels <- prot %>% separate(Lineage, into = levels_vec, sep = ">")
  # Count the occurrance of each group of levels
  protLevels = protLevels %>% group_by_at(levels_vec) %>% summarise(size = n())

  tree <- d3_nest(protLevels, value_cols = "size")

  sb <- sunburst(
    tree,
    #  colors can now also accept a list
    #   specifying range and/or domain of the color scale
    colors = list(
      range=cols,
      domain=legend_items
    ),
    #legendOrder=legend_items, ##### Plotting with legendOrder screws things up
    legend = F,
    width="100%"
  )

  # Add Title to Sunburst
  sb <- htmlwidgets::prependContent(sb, htmltools::h3(query, style = "font-family:Helvetica"))
  sb

}


library(htmltools)
sb_plots <- browsable(
  tagList(
    lapply(
      1:length(queries),
      function(n){
        tags$div(
          style="width:25%; float:left;",
          plt_fun(queries[n])
        )
      }
    )
  )
)


sb_plots

#### Save as png and html ####

# install.packages("webshot")
library(webshot)

setwd("data/figures/")
save_html(sb_plots, "sunburst_queries.html")#, libdir = "data/figures/")

setwd("../../")
webshot("sunburst_queries.html", "data/figures/sunburst_queries.png"
        , zoom = 5) # zoom increases dpi
        #vwidth = 480, vheight = 900)



##### Get Legend #####

# Take lineage column and break into the first to levels
prot <- df %>% select(Lineage)
protLevels <- prot %>% separate(Lineage, into = levels_vec, sep = ">")
# Count the occurrance of each group of levels
protLevels = protLevels %>% group_by_at(levels_vec) %>% summarise(size = n())

tree <- d3_nest(protLevels, value_cols = "size")

sun <- sunburst(
  tree,
  withD3 = T,
  #  colors can now also accept a list
  #   specifying range and/or domain of the color scale
  colors = list(
    range=cols,
    domain=legend_items
  ),
  #legendOrder=legend_items,
  legend = list(w = 150, h = 15, r = 3, s = 3), ##### Plotting with legend screws things up
  width="100%"
)

## Change font to Helvetica
sun2<- htmlwidgets::onRender(
  sun,
  "
    function(el, x) {
    d3.selectAll('.sunburst-legend text').attr('font-family', 'Helvetica');

    // check legend
    d3.select(el).select('.sunburst-togglelegend').property('checked',true);
    // simulate click
    d3.select(el).select('.sunburst-togglelegend').on('click')();
    }

    "
)

save_html(sun2, "sunburst_legend.html")

# Cant zoom in without screwing up font size
webshot("sunburst_legend.html", "data/figures/sunburst_query_legend.png", zoom = 2)
## At this point, I just took a snip of the legend