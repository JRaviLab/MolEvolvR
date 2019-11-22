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
domain_network <- function(prot, column = "DomArch", domains_of_interest, cutoff_type = "Lineage", cutoff = 1, layout = "grid"){
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

  column_name <- as.name(column)
  if(cutoff_type == "Lineage"){
    lin_summary <- prot %>% summ.DA.byLin() %>% summ.DA()
    doms_above_cutoff <- (lin_summary %>% filter(totallin >= cutoff))[[column_name]]
  }
  else if(cutoff_type == "Total Count"){ #Change this type?
    doms_above_cutoff <- (prot %>% total_counts( type ="DA", cutoff = cutoff))[[column_name]]
  }

  prot <- prot[which(prot[[column_name]] %in% doms_above_cutoff),]

  #prot.list=total$arch
  #prot.list <- unlist(lapply(prot.list, function(x) gsub(x,pattern = "\\?",replacement = "X")))
  prot$DomArch.ntwrk <- as_vector(prot %>% select({{column}})) %>% # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  #dom="LTD"  #your domain name here
  domains_of_interest <- c("PspA || PspA_IM30 || PspA\\_IM30")  # your domain name here
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
  if(any(pwise=="")) {
    pwise=pwise[-which((pwise[,1]=="") | pwise[,2]==""),]
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

  #plot.igraph(g,layout = layout_randomly(g),vertex.label.dist=0,asp=0,edge.curved=F, vertex.label.color="black")
  ### Tkplot does not work with shiny
  #tkplot(g, layout=layout_with_gem(g),
  #       vertex.label.dist=1,asp=0,
  #       edge.curved=F, vertex.label.color="black") #simple
}



#pspa_grouped <- pspa %>% group_by(DomArch) %>% summarise(n())
#Counts of all the relations: bind back into df and then can do filtering based off it?
#Else can i