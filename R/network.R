## Function to create a domain network
## Written by: LA & VA
## Modified by: JR, SZC
## Last modified: Sep 24 2019

#################
## Pkgs needed ##
#################
library(tidyverse)
library(igraph)
conflicted::conflict_prefer("filter", "dplyr")

###########################
#### Network FUNCTIONS ####
###########################
domain_network <- function(prot, column = "DomArch", domains_of_interest, cutoff = 1, layout = "grid"){
  #'Domain Network
  #'
  #'This function creates a domain network from the 'DomArch' column.
  #'
  #'A network of domains is returned based on shared domain architectures.
  #'
  #'@param prot A data frame that contains the column 'DomArch'.
  #'@param column Name of column containing Domain architecture from which nodes and edges are generated.
  #'@param cutoff Integer. Only use domains that occur at or above the cutoff
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto"}
  #'@examples domain_network(pspa)

  # by domain networks or all, as required.
  # ye is either all of prot.list or centered on one domain

  column_name <- sym(column)

  prot$DomArch.ntwrk <- as_vector(prot %>% select({{column}})) %>% # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  #dom="LTD"  #your domain name here
  #domains_of_interest <- c("PspA || PspA_IM30 || PspA\\_IM30")  # your domain name here
  #ye=grep(pattern = dom,x = prot.list,value = T)
  #ye=unlist(strsplit(ye,"\\+"))
  ye <- prot %>%
    dplyr::filter(grepl(pattern=domains_of_interest,
                        x=DomArch.ntwrk,
                        ignore.case=T))

  ##Separating column and converting to atomic vector prevents coercion
  ye <- ye$DomArch.ntwrk  %>% str_split(pattern="\\+")

  # Domain network
  domain.list <- ye#strsplit(ye,"\\+") #[[3]] or whatever index of the list(ye) contains domntwk
  if(any(lengths(domain.list)==1)) {
    domain.list=domain.list[-which(lengths(domain.list)==1)]
  }
  te <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
  ye <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n
  pwise <- cbind(te,ye)

  ### Use dom counts from all of prot for node weights
  dom_count <- elements2words(prot,column = {{column}},conversion_type = "da2doms") %>% words2wc()

  ### Use counts from domlist (That filter out the singles) for the edgeweights
  ## Get wordcount for domains (for filtering purposes)
  #dl <- data.frame("Doms" = unlist(domain.list), stringsAsFactors = F)
  #dom_count <- elements2words(dl,column = "Doms",conversion_type = "da2doms") %>% words2wc()


  TotalWordsFreq  = sum(dom_count$freq)
  dom_count <- dom_count %>% mutate("WordsPercentage" = 0) %>% arrange(freq)

  total_counter = 0
  for(x in 1:length(dom_count$freq)){
    total_counter = total_counter + dom_count$freq[x]
    dom_count$WordsPercentage[x] = total_counter/TotalWordsFreq * 100
  }

  doms_above <- dom_count %>% filter(WordsPercentage >= 100 - cutoff)
  if(length(doms_above$words) == 0){
    cutoff_count = 0
  }
  else{
    cutoff_count = doms_above$freq[1]
  }

  # node weights contains the occurances of all domains above the count
  node_weights <- (dom_count %>% filter(freq >= cutoff_count) %>% select(-WordsPercentage))
  ####


  doms_above <- paste(node_weights$words, collapse = "$|^")
  doms_above <- paste0("^",doms_above, "$")

  # Filter out the ones rows containing edges that contain domains below cutoff
  if(any(pwise==""|!grepl(doms_above,pwise))) {
    pwise=pwise[-which((pwise[,1]=="") | pwise[,2]==""
                       | !grepl(doms_above, pwise[,1]) | !grepl(doms_above, pwise[,2]) ),]
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

  #plot.igraph(g,layout = layout_randomly(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black")
  ### Tkplot does not work with shiny
  #tkplot(g, layout=layout_with_gem(g),
  #       vertex.label.dist=1,asp=0,
  #       edge.curved=F, vertex.label.color="black") #simple
}
