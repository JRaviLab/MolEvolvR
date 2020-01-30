# Undirected Graph for GenContext

gc_undirected_network <- function(prot, column = "GenContext", domains_of_interest, cutoff_type = "Lineage", cutoff = 1, layout = "grid"){
  #'Domain Network
  #'
  #'This function creates a domain network from the 'DomArch' column.
  #'
  #'A network of domains is returned based on shared domain architectures.
  #'
  #'@param prot A data frame that contains the column 'DomArch'.
  #'@param column Name of column containing Domain architecture from which nodes and edges are generated.
  #'@param cutoff_type Character. Used to determine how data should be filtered. Either
  #'\itemize{\item "Lineage" to filter domains based off how many lineages the Domain architecture appears in
  #'\item "Total Count" to filter off the total amount of times a domain architecture occurs }
  #'@param cutoff Integer. Only use domains that occur at or above the cutoff for total counts if cutoff_type is "Total Count".
  #'Only use domains that appear in cutoff or greater lineages if cutoff_type is Lineage.
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto"}
  #'@examples domain_network(pspa)

  # by domain networks or all, as required.
  # ye is either all of prot.list or centered on one domain

  column_name <- sym(column)
  if(cutoff_type == "Lineage"){
    lin_summary <- prot %>% summ.DA.byLin() %>% summ.DA()
    doms_above_cutoff <- (lin_summary %>% filter(totallin >= cutoff))[[column]]
  }
  else if(cutoff_type == "Total Count"){ #Change this type?
    GC_above_cutoff <- (prot %>% total_counts( column =  column, cutoff = cutoff))[[column]]
  }

  prot <- prot[which(prot[[as_string(column_name)]] %in% GC_above_cutoff),]

  #prot.list=total$arch
  #prot.list <- unlist(lapply(prot.list, function(x) gsub(x,pattern = "\\?",replacement = "X")))
  prot$GC.ntwrk <- as_vector(prot %>% select({{column}})) %>% # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  ### Do I need a GC of interest?
  # ye <- prot %>%
  #   dplyr::filter(grepl(pattern=domains_of_interest,
  #                       x=DomArch.ntwrk,
  #                       ignore.case=T))

  ## Split the GC column by arrows and '||'
  GC_split <- prot$GC.ntwrk  %>% str_split(pattern="<-|->|\\|\\|")


  #GC.list <- GC_split
  if(any(lengths(GC_split)==1)) {
    GC_split = GC_split[-which(lengths(GC_split)==1)]
  }
  lead <- unlist(lapply(GC_split, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
  follow <- unlist(lapply(GC_split, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n
  pwise <- cbind(lead,follow)
  if(any(pwise==""|pwise=="?")) {
    pwise=pwise[-which((pwise[,1]=="") | (pwise[,1]=="?") | (pwise[,2]=="")| (pwise[,2]=="?")),]
  }
  pwise.ass <- sapply(1:length(pwise[,1]), function(x) paste(pwise[x,],collapse = "->"))
  e.sz <- sort(table(pwise.ass),decreasing = T)  # Edge weights
  v.sz <- sort(table(pwise),decreasing = T)   # Vertex weights
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
