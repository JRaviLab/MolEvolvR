## Function to create a domain network
## Written by: LA & VA
## Modified by: JR, SZC
## Last modified: Sep 24 2019

#################
## Pkgs needed ##
#################
library(igraph)
library(tidyverse)

###########################
#### CLEANUP FUNCTIONS ####
###########################
#'Cleanup Species
#'
#'This function creates a domain network from the 'DomArch' column.
#'
#'A network of domains is returned based on shared domain architectures.
#'
#'@param prot A data frame that contains the column 'DomArch'.
#'@examples domain_network(pspa)
domain_network <- function(prot){
# by domain networks or all, as required.
# ye is either all of prot.list or centered on one domain

prot.list <- prot$DomArch # testing with prot$DomArch.orig
prot.list <- unlist(lapply(prot.list, function(x) gsub(x,pattern = "\\?",replacement = "X")))
domains_of_interest <- c("PspA || PspA_IM30 || PspA\\_IM30")  # your domain name here
ye <- prot.list %>%
  grep(pattern=domains_of_interest, value = T)
ye <- unlist(strsplit(ye,"\\+"))

# Domain network
domain.list <- strsplit(ye,"\\+")
if(any(lengths(domain.list)==1)) {
  domain.list=domain.list[-which(lengths(domain.list)==1)]
  }
te <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y])))
ye <- unlist(lapply(domain.list, function(x) sapply(1:(length(x)-1), function(y) x[y+1])))
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
V(g)$color <- rainbow(5)[round((V(g)$size-min(V(g)$size))/(max(V(g)$size)-min(V(g)$size))*4+1)]
V(g)$frame.color <- V(g)$color
# scaling edge width
E(g)$width <- e.sz
E(g)$width <- ifelse(log(E(g)$width)==0, .3,log(E(g)$width))
# coloring edges by width
ew <- c(2.7,4.5)
E(g)$color <- sapply(E(g)$width,
                  function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color="cadetblue" else if(x>ew[2]) E(g)$color="maroon" else E(g)$color="gray55")

tkplot(g, layout=layout_with_gem(g),
       vertex.label.dist=1,asp=0,
       edge.curved=F, vertex.label.color="black") #simple
}