source("plotting.R")

par(mfrow=c(4,2))

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

#lin_counts <- prot_all %>% group_by(Lineage) %>% summarise(count = n()) %>% arrange(-count)

# # Take lineage column and break into the first to levels
# prot_all <- prot_all %>% select(Lineage)
# protLevels <- prot_all %>% separate(Lineage, into = c("Level1","Level2"), sep = ">")
# # Count the occurrance of each group of levels
# protLevels = protLevels %>% group_by_at(levels_vec) %>% summarise(size = n())
#
# tree <- d3_nest(protLevels, value_cols = "size")


df <- prot_all

rownames(df) <- NULL

legend_items <- unique(unlist(strsplit(df$Lineage, ">")))

cols <- sample(colorRampPalette(brewer.pal(12, 'Set3'))(length(legend_items)))

queries = c("PspA", "Snf7", "Psp-AA", "PspB","PspC",  "PspN", "LiaI-LiaF-TM", "Toast-rack")

# plot each lineage
plt_fun <- function(query){
  # filter for dt
  plt_df <- df %>%
    filter(grepl(query, DomArch, ignore.case = T)) %>%
    select(-DomArch)

  # Take lineage column and break into the first to levels
  prot <- plt_df %>% select(Lineage)
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
    legendOrder=legend_items,
    width="100%"
  )
}


library(htmltools)
browsable(
  tagList(
    lapply(
      1:2,
      function(n){
        tags$div(
          style="width:48%; float:left;",
          plt_fun(queries[n])
        )
      }
    )
  )
)



all_sun <- lineage_sunburst(all)

pspa_sun <- lineage_sunburst(pspa, "Lineage", "sunburst", 2, legendOrder = legend_order)

pspb_sun <- lineage_sunburst(pspb, "Lineage", "sunburst", 2)
