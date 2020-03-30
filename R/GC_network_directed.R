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
                                 cutoff = 40,
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
  #'@examples gc_directed_network(pspa, column = "GenContex", cutoff = 55)

  column_name <- sym(column)


  # Identify cutoff
  ### Use GC counts from all of prot for node weights
  GC_count <- elements2words(prot,column = {{column}},conversion_type = "gc2da") %>% words2wc()

  # Get Cummulative Percentage of Each GC element
  TotalWordsFreq  = sum(GC_count$freq)
  GC_count <- GC_count %>% mutate("WordsPercentage" = 0) %>% arrange(freq)

  total_counter = 0
  for(x in 1:length(GC_count$freq)){
    total_counter = total_counter + GC_count$freq[x]
    GC_count$WordsPercentage[x] = total_counter/TotalWordsFreq * 100
  }

  GC_above <- GC_count %>% filter(WordsPercentage >= 100 - cutoff)
  if(length(GC_above$words) == 0){
    cutoff_count = 0
  }
  else{
    cutoff_count = GC_above$freq[1]
  }

  # node weights contains the occurances of all domains above the count
  node_weights <- (GC_count %>% filter(freq >= cutoff_count) %>% select(-WordsPercentage))
  ####


  GC_above <- paste(node_weights$words, collapse = "$|^")
  GC_above <- paste0("^",GC_above, "$")


  #### rev_gc(gc_row)
  #### Reverse the order of the split up GC elements whenever a change in transcription direction is found?
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

  # Prime for splitting
  gencontext <- gsub(pattern = ">",replacement = ">|",x = gencontext)
  gencontext <- gsub(pattern = "<",replacement = "|<",x = gencontext)

  # Split GC
  gen_split <- gencontext %>% strsplit("\\|")

  # Remove non edges
  if(any(lengths(gen_split)==1)) {
    gen_split=gen_split[-which(lengths(gen_split)==1)]
  }

  # Apply rev_gc to all rows
  gen_split_reoriented <- map(.x = gen_split, .f = rev_gc)

  # Remove all arrows
  gen_split_reoriented <- map2(.x = gen_split_reoriented,.y = "<-|->", .f = str_remove_all)

  te <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
  ye <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n

  # Edge list
  pwise <- cbind(te,ye)

  # Remove rows that are not in cutoff or have bad characters
  if(any(pwise==""|pwise=="?"| !grepl(GC_above,pwise))) {
    pwise=pwise[-which((pwise[,1]=="") | (pwise[,1]=="?") | (pwise[,2]=="")| (pwise[,2]=="?")
                       | !grepl(GC_above, pwise[,1]) | !grepl(GC_above, pwise[,2]) ),]
  }



  # First create list graph from edge list, then add the
  # vertices that are above the cutoff but not in the graph

  if(length(pwise) != 0){
    # Original
    pwise.ass <- sapply(1:length(pwise[,1]), function(x) paste(pwise[x,],collapse = "->"))
    e.sz <- sort(table(pwise.ass),decreasing = T)   # Counts of each edge (edge weight)
    v.sz <- sort(table(pwise),decreasing = T)       # Counts of each vertex in the edges -- May be obsolete if i just r
    # Relace with al the node weights
    pwise <- strsplit(names(e.sz), "\\->")
    pwise <- cbind(unlist(lapply(pwise, function(x) x[1])),unlist(lapply(pwise, function(x) x[2])))

    g <- graph_from_edgelist(pwise, directed=TRUE)

    ### Handle Edge Weights and Color
    # scaling edge width
    E(g)$width <- e.sz
    E(g)$width <- ifelse(log(E(g)$width)==0, .3,log(E(g)$width))
    # coloring edges by width
    ew <- c(2.7,4.5)
    E(g)$color <- sapply(E(g)$width,
                         function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color="cadetblue" else if(x>ew[2]) E(g)$color="maroon" else E(g)$color="gray55")
  }
  else{
    g <- make_empty_graph()
  }

  # Go through all vertices above cutoff, add them to graph if they don't exist
  # Assign them their respective weights with V(g)$size <- v.sz[V(g)$name]

  node_weights <- pivot_wider(node_weights, names_from = words, values_from = freq)

  # vertices to be added -- No edges in graph
  v_to_add <- colnames(node_weights)[which(!(colnames(node_weights) %in% V(g)$name))]

  # Add isolated vertices
  g <- g + vertices(v_to_add)

  # Assign Weights to each vertex
  V(g)$size <- as.integer(node_weights[V(g)$name])

  # if only one node, set to fixed size Or if min and max are same number
  if(length(V(g)$size) == 1 | min(V(g)$size) == max(V(g)$size)  )
  {
    V(g)$size <- 25
    V(g)$color <- rainbow(1, alpha = .5)
    V(g)$frame.color <- V(g)$color
  }

  # Resize based on difference between sizes of max and size of min
  # (length(V(g)$size) >1)
  else{
    ## Below scaling will not work if only one vertex
    V(g)$size <- (V(g)$size-min( V(g)$size) )/(max(V(g)$size)-min(V(g)$size))*20+10 # scaled by degree

    v_size <-(V(g)$size)
    # setting vertex color by size
    V(g)$color <- rainbow(5,alpha = .5)[round( (v_size-min(v_size))/(max(v_size)-min(v_size)*4+1))]
    V(g)$frame.color <- V(g)$color
  }

  switch(layout,
         "random" = plot.igraph(g,layout = layout_randomly(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "grid" = plot.igraph(g,layout = layout_on_grid(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "circle" = plot.igraph(g,layout = layout.circle(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black"),
         "auto" =plot.igraph(g,layout = layout.auto(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black")
  )
}
