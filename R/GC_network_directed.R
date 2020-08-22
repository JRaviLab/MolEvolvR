library(tidyverse)
library(rlang)
library(igraph)
library(visNetwork)
conflicted::conflict_prefer("filter", "dplyr")

GenContextNetwork <- function(prot,domains_of_interest ,column = "GenContext",
                              cutoff = 40,
                              layout = "grid",
                              directed = TRUE){
  #'Genomic Context Directed Network
  #'
  #'This function creates a Genomic Context network from the 'GenContext' column.
  #'
  #'A network of Genomic Context is returned.
  #'
  #'@param prot A data frame that contains the column 'GenContext'.
  #'@param column Name of column containing Genomic Context from which nodes and edges are generated.
  #'@param domains_of_interest Character vector of domains of interest.
  #'@param cutoff Integer. Only use GenContexts that occur at or above the cutoff percentage for total count
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto" \item "nice"}
  #'@param directed Is the network directed?
  #'@examples gc_directed_network(pspa, column = "GenContex", cutoff = 55)

  column_name <- sym(column)


  # Perform cutoff through total_counts
  prot_tc <- prot %>% total_counts(column =  column, cutoff = cutoff)

  within_list <- prot_tc %>% select({{column_name}}) %>% distinct()
  within_list <- pull(within_list, {{column_name}})

  prot <- prot %>% filter({{column_name}} %in% within_list)

  # Below is obsolete kinda
  prot$DomArch.ntwrk <- as_vector(prot %>% select({{column}})) %>%
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  # # Make sure all observations contain the domain of interest
  # domains_of_interest_regex = paste(domains_of_interest, collapse = "|")
  # ye <- prot %>%
  #   dplyr::filter(grepl(pattern=domains_of_interest_regex,
  #                       x=DomArch.ntwrk,
  #                       ignore.case=T))

  # ##Separating column and converting to atomic vector prevents coercion
  # ye <- ye$DomArch.ntwrk  %>% str_split(pattern="\\+")
  #
  # GC.list <- ye
  # if(any(lengths(GC.list)==1)) {
  #   GC.list=GC.list[-which(lengths(GC.list)==1)]
  # }


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

  # Get domain counts before eliminating domarchs with no edges
  wc = elements2words(prot = prot, column =  column, conversion_type = "gc2da") %>% words2wc()
  nodes = data.frame(id = wc$words, label = wc$words, size = wc$freq)

  max_size = max(nodes$size)
  min_size = min(nodes$size)
  nodes <- nodes %>% mutate(size = (size-min_size)/((max_size-min_size)) *20+10)
  max_font_size = 43
  nodes <- nodes %>% mutate(font.size = purrr::map(size, function(x) min(x*2 , max_font_size)))

  max_size = max(nodes$size)
  min_size = min(nodes$size)
  nodes <- nodes %>% mutate( color = unlist(purrr::map(nodes$size ,
                                             function(x)
                                             {
                                               return(rainbow(5,alpha = .5)[round( (x-min_size)/(max_size-min_size)*4+1 )])
                                             }
  ))  )
  nodes$frame.color <- nodes$color



  gencontext <- prot$DomArch.ntwrk

  # Split GC
  gencontext <- gsub(pattern = ">",replacement = ">|",x = gencontext)
  gencontext <- gsub(pattern = "<",replacement = "|<",x = gencontext)

  gen_split <- gencontext %>% strsplit("\\|")

  # Remove non edges
  if(any(lengths(gen_split)==1)) {
    gen_split=gen_split[-which(lengths(gen_split)==1)]
  }

  # Apply rev_gc to all rows
  gen_split_reoriented <- map(.x = gen_split, .f = rev_gc)

  # Remove all arrows
  gen_split_reoriented <- map2(.x = gen_split_reoriented,.y = "<-|->", .f = str_remove_all)

  from <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y]))) #list elements 1 through n-1
  to <- unlist(lapply(gen_split_reoriented, function(x) sapply(1:(length(x)-1), function(y) x[y+1]))) #list elements 2 through n

  # Edge list
  pwise <- cbind(from,to)

  # Remove rows that have bad characters
  if(any(pwise==""|pwise=="X"| pwise== "X(s)")) {
    pwise=pwise[-which((pwise[,1]=="") | (pwise[,1]=="X")|(pwise[,1]=="X(s)")  | (pwise[,2]=="")| (pwise[,2]=="X") |(pwise[,2]=="X(s)")  ),]
  }

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

  if(directed){
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px')  %>%
      visEdges(arrows = 'to', smooth =T)
  }
  else
  {
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px')  %>%
      visEdges( smooth =T)
  }

  vg <- vg %>%
    visOptions(highlightNearest = TRUE) %>%
    visLegend()

  vg <- switch(layout,
               "nice" = visIgraphLayout(vg, "layout_nicely" ),
               "random" = visIgraphLayout(vg, "layout_randomly"),
               "grid" = visIgraphLayout(vg, "layout_on_grid"),
               "circle" = visIgraphLayout(vg, "layout.circle"),
               "auto" =  visIgraphLayout(vg, "layout.auto")
  )
  vg


}
