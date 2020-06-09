source("plotting.R")

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
cols <- sample(colorRampPalette(brewer.pal(9, 'Set1'))(length(legend_items)))
# cols <- sample(colorRampPalette(brewer.pal(12, 'Paired'))(length(legend_items)))
# cols <- sample(colorRampPalette(brewer.pal(12, 'Pastel2'))(length(legend_items)))

queries = c("PspA", "Snf7", "PspB","PspC", "PspM",  "PspN", "LiaI-LiaF-TM", "Toast-rack")

# df$Lineage <- str_replace_all(df$Lineage, "-", "_")
#
# df$relationship <- str_replace_all(df$Lineage, ">", "-")

# plot each lineage
plt_fun <- function(query){
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

save_html(sb_plots, "sunburst_queries.html")

webshot("sunburst_queries.html", "sunburst_queries.png")
        #vwidth = 480, vheight = 900)



##### Get Legend #####

# Take lineage column and break into the first to levels
prot <- df %>% select(Lineage)
protLevels <- prot %>% separate(Lineage, into = levels_vec, sep = ">")
# Count the occurrance of each group of levels
protLevels = protLevels %>% group_by_at(levels_vec) %>% summarise(size = n())

tree <- d3_nest(protLevels, value_cols = "size")

sunburst(
  tree,
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

## At this point, I just took a snip of the legend
