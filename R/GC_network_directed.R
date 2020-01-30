#################
## Pkgs needed ##
#################
library(tidyverse)
library(igraph)
conflicted::conflict_prefer("filter", "dplyr")

###########################
#### Network FUNCTIONS ####
###########################
gc_directed_network <- function(prot, column = "GenContext",
                       cutoff_type = "Lineage", cutoff = 40,
                       layout = "grid"){
  #'Genomic Context Directed Network
  #'
  #'This function creates a Genomic Context network from the 'GenContext' column.
  #'
  #'A network of Genomic Context is returned.
  #'
  #'@param prot A data frame that contains the column 'GenContext'.
  #'@param column Name of column containing Genomic Context from which nodes and edges are generated.
  #'@param cutoff_type Character. Used to determine how data should be filtered. Either
  #'\itemize{\item "Lineage" to filter domains based off how many lineages the Genomic Context appears in
  #'\item "Total Count" to filter off the total amount of times a domain architecture occurs }
  #'@param cutoff Integer. Only use GenContexts that occur at or above the cutoff percentage for total counts if cutoff_type is "Total Count".
  #'Only use domains that appear in cutoff or greater lineages if cutoff_type is Lineage.
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto"}
  #'@examples gc_directed_network(pspa, cutoff_type = "Total Count", cutoff = 55)

  # by domain networks or all, as required.
  # ye is either all of prot.list or centered on one domain

  column_name <- sym(column)
  if(cutoff_type == "Lineage"){
    # Cutoff is determined by how many different lineages the GC appears in
    lin_summary <- prot %>% summ.GC.byLin() %>% group_by(GenContext) %>% summarize(totallin = n())
    GC_above_cutoff <- (lin_summary %>% filter(totallin >= cutoff))[[column]]
  }
  else if(cutoff_type == "Total Count"){
    # Cutoff is determined by the overall amount of times a GC appears
    GC_above_cutoff <- (prot %>% total_counts( column =  column, cutoff = cutoff))[[column]]
  }

  prot <- prot[which(prot[[{{column_name}}]] %in% GC_above_cutoff),]

  rev_gc <- function(x){
    direction <- 0
    op_dir <- "<-|->"
    # if direction is 1, '->' are the correct arrows
    # if direction is -1, '<-'
    last <- ""
    for(i in 1:length(x)){
      if(direction == 0 && grepl("<-|->",x[i])){
        # Set initial direction of the GC
        if(grepl("^<-", x[i])){
          op_dir <- "->"
          direction = -1
        }
        else{
          op_dir <- "<-"
          direction = 1
        }
      }
      else if(grepl(op_dir, x[i])){
        # Arrow going in the opposite direction has been found
        # Reverse order of the remaining list--- Do I need to add a " " before this?
        x[i:length(x)] = rev(x[i:length(x)])
        if(direction == 1){
          op_dir <- "->"
          direction = -1
        }
        else if(direction == -1){
          op_dir <- "<-"
          direction = 1
        }
      }
    }
    return(x)
  }

  gencontext <- prot$GenContext

  gencontext <- gsub(pattern = ">",replacement = ">|",x = gencontext)

  gencontext <- gsub(pattern = "<",replacement = "|<",x = gencontext)


  gen_split <- gencontext %>% strsplit("\\|")

  if(any(lengths(gen_split)==1)) {
    gen_split=gen_split[-which(lengths(gen_split)==1)]
  }

  gen_split_reoriented <- map(.x = gen_split, .f = rev_gc)


  gen_split_reoriented <- map2(.x = gen_split_reoriented,.y = "<-|->", .f = str_remove_all)


  te <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
  ye <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n
  pwise <- cbind(te,ye)
  if(any(pwise==""|pwise=="?")) {
    pwise=pwise[-which((pwise[,1]=="") | (pwise[,1]=="?") | (pwise[,2]=="")| (pwise[,2]=="?")),]
  }



  pwise.ass <- sapply(1:length(pwise[,1]), function(x) paste(pwise[x,],collapse = "->"))
  e.sz <- sort(table(pwise.ass),decreasing = T)
  v.sz <- sort(table(pwise),decreasing = T)
  pwise <- strsplit(names(e.sz), "\\->")
  pwise <- cbind(unlist(lapply(pwise, function(x) x[1])),unlist(lapply(pwise, function(x) x[2])))

  # Plotting the network
  g <- graph_from_edgelist(pwise, directed=TRUE)
  # scaling vertex size
  V(g)$size <- v.sz[V(g)$name]
  V(g)$size <- (V(g)$size-min(V(g)$size))/(max(V(g)$size)-min(V(g)$size))*20+10 # scaled by degree
  # setting vertex color by size
  V(g)$color <- rainbow(5,alpha = .5)[round((V(g)$size-min(V(g)$size))/(max(V(g)$size)-min(V(g)$size))*4+1)]
  V(g)$frame.color <- V(g)$color
  # scaling edge width
  E(g)$width <- e.sz
  E(g)$width <- ifelse(log(E(g)$width)==0, .3,log(E(g)$width))
  # coloring edges by width
  ew <- c(2.7,4.5)
  E(g)$color <- sapply(E(g)$width,
                       function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color="cadetblue" else if(x>ew[2]) E(g)$color="maroon" else E(g)$color="gray55")


  switch(layout,
         "random" = plot.igraph(g,layout = layout_randomly(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "grid" = plot.igraph(g,layout = layout_on_grid(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "circle" = plot.igraph(g,layout = layout.circle(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "auto" =plot.igraph(g,layout = layout.auto(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black")
  )
}
