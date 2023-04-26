#by domain networks or all, as required.  ye is either all of prot.list or centered on one domain
prot.list=total$arch
prot.list=unlist(lapply(prot.list, function(x) gsub(x,pattern = "\\?",replacement = "X")))
dom="LTD"  #your domain name here
ye=grep(pattern = dom,x = prot.list,value = T)
ye=unlist(strsplit(ye,"\\+"))
 
#domain network
dom.list=strsplit(ye,"\\+")
if(any(lengths(dom.list)==1)) dom.list=dom.list[-which(lengths(dom.list)==1)]
te=unlist(lapply(dom.list, function(x) sapply(1:(length(x)-1), function(y) x[y])))
ye=unlist(lapply(dom.list, function(x) sapply(1:(length(x)-1), function(y) x[y+1])))
pwise=cbind(te,ye)
if(any(pwise=="")) pwise=pwise[-which((pwise[,1]=="") | pwise[,2]==""),]
pwise.ass=sapply(1:length(pwise[,1]), function(x) paste(pwise[x,],collapse = "->"))
e.sz=sort(table(pwise.ass),decreasing = T)
v.sz=sort(table(pwise),decreasing = T)
pwise=strsplit(names(e.sz), "\\->")
pwise=cbind(unlist(lapply(pwise, function(x) x[1])),unlist(lapply(pwise, function(x) x[2])))
 
suppressPackageStartupMessages(library("igraph"))
g = graph_from_edgelist(pwise, directed=TRUE)
#scaling vertex size
V(g)$size=v.sz[V(g)$name]
V(g)$size=(V(g)$size-min(V(g)$size))/(max(V(g)$size)-min(V(g)$size))*20+10 # scaled by degree
#setting vertex color by size
V(g)$color= rainbow(5)[round((V(g)$size-min(V(g)$size))/(max(V(g)$size)-min(V(g)$size))*4+1)]
V(g)$frame.color=V(g)$color
#scaling edge width
E(g)$width=e.sz
E(g)$width=ifelse(log(E(g)$width)==0, .3,log(E(g)$width))
#coloring edges by width
ew=c(2.7,4.5)
E(g)$color=sapply(E(g)$width,
                  function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color="cadetblue" else if(x>ew[2]) E(g)$color="maroon" else E(g)$color="gray55")
 
tkplot(g, layout=layout_with_gem(g), vertex.label.dist=1,asp=0,edge.curved=F, vertex.label.color="black") #simple