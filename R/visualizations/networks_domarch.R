## Function to create a domain network
## Written by: LA & VA
## Modified by: JR, SZC
## Last modified: Sep 24 2019
#################
## Pkgs needed ##
#################
library(tidyverse)
library(rlang)
library(igraph)
library(visNetwork)
conflicted::conflict_prefer("filter", "dplyr")
###########################
#### Network FUNCTIONS ####
###########################

##### !!!!! CHANGE TO MAKE OUR OWN !!!!!!

domain_network <- function(prot, column = "DomArch", domains_of_interest, cutoff = 70, layout = "nice", query_color = adjustcolor("green", alpha.f = .5)){
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
  print(domains_of_interest)
  tryCatch({
  column_name <- sym(column)

  prot_tc <- prot %>% total_counts(column =  column, cutoff = cutoff, RowsCutoff = F, digits = 5)

  within_list <- prot_tc %>% select({{column_name}}) %>% distinct()
  within_list <- pull(within_list, {{column_name}})

  # nvm, it's making sure that only Domains that are in the tc cutoff range are kept
  prot <- prot %>% filter({{column_name}} %in% within_list)

  ####### Below should be part of the standardized cleanup process
  prot$DomArch.ntwrk <- as_vector(prot %>% select({{column}})) %>% # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")


  #dom="LTD"  #your domain name here
  #domains_of_interest <- c("PspA || PspA_IM30 || PspA\\_IM30")  # your domain name here
  #ye=grep(pattern = dom,x = prot.list,value = T)
  #ye=unlist(strsplit(ye,"\\+"))
  domains_of_interest_regex = paste(domains_of_interest, collapse = "|")
  domain.list <- prot %>%
    dplyr::filter(grepl(pattern=domains_of_interest_regex,
                        x=DomArch.ntwrk,
                        ignore.case=T, perl = T))

 
  ##Separating column and converting to atomic vector prevents coercion
  domain.list <- domain.list$DomArch.ntwrk  %>% str_split(pattern="\\+")
  # Get domain counts before eliminating domarchs with no edges
  wc = elements2words(prot = prot, column =  column, conversion_type = "da2doms") %>% words2wc()
  wc = pivot_wider(wc, names_from = words, values_from = freq)

  # Remove all isolated domarchs, such that an adjacency list can easily be constructed
  singletons <- domain.list[which(lengths(domain.list)==1)] %>% unique()
  domain.list = domain.list[-which(lengths(domain.list)==1)]
  # This is where we know if the adjacency list is empty
  if(length(domain.list) == 0)
  {
    g <- make_empty_graph()
    # Add nodes included in domains of interest
    g <- g + vertices(singletons)
    V(g)$size <- length(singletons)

    if(length(V(g)$size) == 1 || min(V(g)$size) == max(V(g)$size)  )
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



  }

  else
  {
    te <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
    ye <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n
    pwise <- cbind(te,ye)
    if(any(pwise=="")) {
      pwise=pwise[-which((pwise[,1]=="") | pwise[,2]==""),]
    }
    pwise.ass <- sapply(1:length(pwise[,1]), function(x) paste(pwise[x,],collapse = "->"))
    e.sz <- sort(table(pwise.ass),decreasing = T)
    # v.sz <- sort(table(pwise),decreasing = T)
    e.sz <- sort(table(pwise.ass),decreasing = T)  # Edge weights
    # v.sz <- sort(table(pwise),decreasing = T)   # Vertex weights
    pwise <- strsplit(names(e.sz), "\\->")
    pwise <- cbind(unlist(lapply(pwise, function(x) x[1])),unlist(lapply(pwise, function(x) x[2])))

    # Plotting the network
    g <- graph_from_edgelist(pwise, directed=TRUE)
    # scaling vertex size
    #V(g)$size <- v.sz[V(g)$name]


    # # Add query domains if not already present b/c they have not adjacency add them to the graph if not already present
    for(domain in singletons)
    {
      if (! domain %in% V(g)$name){
        g <- g + vertex(domain)
      }
    }

    # Make sure X does not appear
    # if(g)
    if("X" %in% V(g)$name){
      g = delete_vertices(g, "X")
    }
    V(g)$size <- as.numeric(wc[V(g)$name])

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
                         function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color=adjustcolor("cadetblue", alpha.f = .7) else if(x>ew[2]) E(g)$color=adjustcolor("maroon", alpha.f = .5) else E(g)$color="gray55")
  }


  #   # V(g)[which(V(g)$name %in% domains_of_interest)]$size <- as.numeric(wc[domains_of_interest])
  V(g)[which(V(g)$name %in% domains_of_interest)]$color = query_color

  vis_g = visIgraph(g, type ="full")

  max_font_size = 43

  vis_g$x$nodes$font.size <- purrr::map(vis_g$x$nodes$size, function(x) min(x*2 , max_font_size))

  vis_g <- switch(layout,
                  "nice" = visIgraphLayout(vis_g,"layout_nicely", ),
                  "random" = visIgraphLayout(vis_g, "layout_randomly"),
                  "grid" = visIgraphLayout(vis_g, "layout_on_grid"),
                  "circle" = visIgraphLayout(vis_g, "layout.circle"),
                  "auto" =  visIgraphLayout(vis_g, "layout.auto")
  )
  vis_g <- vis_g %>%
    visOptions(highlightNearest = TRUE)
  },
  error = {
    vis_g <- "error"
  },
  finally = {
    return(vis_g)
  })
                                        
}

BinaryDomainNetwork <- function(prot, column = "DomArch", domains_of_interest, cutoff = 70,
                                layout = "nice", query_color = adjustcolor("yellow", alpha.f = .5),
                                partner_color = adjustcolor("skyblue", alpha.f = .5),
                                border_color = adjustcolor("grey", alpha.f = .8),
                                IsDirected = T
){
  #'Domain Network
  #'
  #'This function creates a domain network from the 'DomArch' column.
  #'
  #'Domains that are part of the 'domains_of_interest' are a different node color than the other domains.
  #'
  #'A network of domains is returned based on shared domain architectures.
  #'
  #'@param prot A data frame that contains the column 'DomArch'.
  #'@param column Name of column containing Domain architecture from which nodes and edges are generated.
  #'@param cutoff Integer. Only use domains that occur at or above the cutoff for total counts if cutoff_type is "Total Count".
  #'Only use domains that appear in cutoff or greater lineages if cutoff_type is Lineage.
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto"}
  #'@param query_color Color that the nodes of the domains in the domains_of_interest vector are colored
  #'@param partnercolor Color that the nodes that are not part of the domains_of_interest vector are colored
  #'@param IsDirected Is the network directed? Set to false to eliminate arrows
  #'@examples domain_network(pspa)
  # by domain networks or all, as required.
  print(domains_of_interest)

  column_name <- sym(column)

  prot_tc <- prot %>% total_counts(column =  column, cutoff = cutoff, RowsCutoff = F, digits = 5)

  within_list <- prot_tc %>% select({{column_name}}) %>% distinct()
  within_list <- pull(within_list, {{column_name}})

  prot <- prot %>% filter({{column_name}} %in% within_list)

  ####### Below should be part of the standardized cleanup process
  prot$DomArch.ntwrk <- as_vector(prot %>% select({{column}})) %>% # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  domains_of_interest_regex = paste(domains_of_interest, collapse = "|")
  # domain.list <- prot %>%
  #   dplyr::filter(grepl(pattern=domains_of_interest_regex,
  #                       x=DomArch.ntwrk,
  #                       ignore.case=T))
  domain.list = prot

  ##Separating column and converting to atomic vector prevents coercion
  domain.list <- domain.list$DomArch.ntwrk  %>% str_split(pattern="\\+")

  # Get domain counts before eliminating domarchs with no edges
  wc = elements2words(prot = prot, column =  column, conversion_type = "da2doms") %>% words2wc()

  nodes = data.frame(id = wc$words, label = wc$words, size = wc$freq) %>%
    mutate(group = purrr::map(id,
                              function(x)
                                ifelse(x %in% domains_of_interest, "Query", "Partner")))

  max_size = max(nodes$size)
  min_size = min(nodes$size)
  nodes <- nodes %>% mutate(size = (size-min_size)/((max_size-min_size)) *20+10)
  max_font_size = 43
  nodes <- nodes %>% mutate(font.size = purrr::map(size, function(x) min(x*2 , max_font_size)))

  domain.list = domain.list[-which(lengths(domain.list)==1)]

  if(length(domain.list) != 0)
  {
    from <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
    to <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n
    pwise <- cbind(from,to)
    if(any(pwise=="")) {
      pwise=pwise[-which( (pwise[,1]=="") | (pwise[,2]=="") ),]
    }
    if(any(pwise == "X"))
      pwise=pwise[-which( (pwise[,1]=="X") | (pwise[,2]=="X") ),]

    edges <- data.frame(from = pwise[,1], to = pwise[,2]) %>%
      group_by(from, to) %>% summarize(width = n())
    edges <- edges %>% mutate(width = ifelse(width==1, .3, log(width)))
    ew <- c(2.7,4.5)

    ColorEdges <- function(x)
    {
      if(x>=ew[1] && x<=ew[2])
      {
        adjustcolor("cadetblue", alpha.f = .7)
      }
      else if(x>ew[2])
      {
        adjustcolor("maroon", alpha.f = .5)
      }
      else
      {
        "gray55"
      }
    }

    edges <- edges %>% mutate(color = unlist(purrr::map(width, ColorEdges)) )

  }

  if(IsDirected){
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px')  %>%
      visEdges(arrows = 'to', smooth =T)
  }
  else
  {
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px')  %>%
      visEdges(smooth =T)
  }
  vg <- vg %>%   visGroups(groupname = "Query", color = query_color) %>%
    visGroups(groupname = "Partner", color = partner_color) %>%
    visOptions(highlightNearest = TRUE) %>%
      visLegend(position = "right",width = .1)

  vg <- switch(layout,
               "nice" = visIgraphLayout(vg, "layout_nicely" ),
               "random" = visIgraphLayout(vg, "layout_randomly"),
               "grid" = visIgraphLayout(vg, "layout_on_grid"),
               "circle" = visIgraphLayout(vg, "layout.circle"),
               "auto" =  visIgraphLayout(vg, "layout.auto")
  )
  vg                                               
}


