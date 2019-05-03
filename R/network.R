## Last modified: April 08, 2019

#############
# basic     #
#############

# install.packages("igraph")
# library(igraph); library(rgl)

te=data.frame(name=c("A","B","C","D","E"),
              sz=c(48,33,45,34,21),
              attr=c("B","S","B","S","B")
              )
ye=data.frame(from=c("B","C","C","D","D","E"),
              to=c("A","B","A","A","B","A"),
              e1=c(4,5,5,2,1,1),
              e2=c(4,5,5,4,2,3))
g = graph_from_data_frame(ye, vertices = te, directed=TRUE)

# simple graph
plot(g, vertex.size=te$sz, edge.width=ye$e1,
     vertex.color=rainbow(length(te$name)),edge.color=rainbow(8)[ye$e2])
# interactive graph
tkplot(g, vertex.size=te$sz, edge.width=ye$e1,
     vertex.color=rainbow(length(te$name)),edge.color=rainbow(8)[ye$e2])
#3D
rglplot(g, vertex.size=te$sz, edge.width=ye$e1,
        vertex.color=rainbow(length(te$name)),edge.color=rainbow(8)[ye$e2])
#Kamada Kawai
plot(g, vertex.size=te$sz, edge.width=ye$e1,
     vertex.color=rainbow(length(te$name)),edge.color=rainbow(8)[ye$e2],
     layout=layout_with_kk)

#random scale-free graph
g=barabasi.game(100)

n=15 # first circle
m=10 # second circle
t=seq(0,2*pi-2*pi/n,length=n)
x=rep(10*cos(t), m)
y=rep(10*sin(t), m)

t=seq(0,2*pi-2*pi/m,length=m)

z=40*cos(t)
z=as.vector(sapply(z,function(x) rep(x,n)))
x=x+z
z=40*sin(t)
z=as.vector(sapply(z,function(x) rep(x,n)))
y=y+z

V(g)$x=x
V(g)$y=y

plot(g, vertex.label.dist=0, vertex.size=3, edge.curved=F, asp=1)

############
#reading domain graph
graph_d<- read.csv("graph-for-network")
#clean up
graph_d=data.frame(graph_d)
#converting to format
colnames(graph_d)=c("V1","V2","V3")

graph_d$V1=sapply(graph_d$V1,as.character)
graph_d$V2=sapply(graph_d$V2,as.character)
graph_d$V3=sapply(graph_d$V3,as.character)
graph_d$V3=sapply(graph_d$V3,as.numeric)

graph_d=graph_d[-c(1,2,3,4),]

write.csv(graph_d, "graph-for-network",row.names = F) #writing back cleaned up file

#reading nodes
nodes = read.csv("nodes")
nodes = data.frame(nodes)
nodes$V1=sapply(nodes$V1,as.character)
nodes$V2=sapply(nodes$V2,as.character)
nodes$V2=sapply(nodes$V2,as.numeric)

write.csv(nodes, "nodes",row.names = F) #writing cleaned up file

# making preliminary graph
g = graph_from_data_frame(graph_d, directed=TRUE, vertices = nodes)

#checking out graph
V(g) #vertices

#checking out degree of some node AIRS
te=degree(g)
te["AIRS"]
te = sort(te) # sort by
#check nodes by degree
V(g)[degree(g)==1]

deg=degree(g) # degree of vertices
deg=hist(deg, breaks =seq(0,max(deg),1)) #histogram of degree
plot(seq(1,length(deg$counts),1),deg$counts, log="xy", pch=16) #log scale
plot(seq(1,length(deg$counts),1),deg$counts, pch=16) #plain

fit_power_law(degree(g)+1) # checking power law stats



diameter(g,directed = T) #diameter of graph

distances(g, v="MPTase",to=V(g)) #checking distances of all nodes from MPTase
hist(distances(g, v="MuF",to=V(g)),breaks = c(-.5,.5, 1.5,2.5,3.5,4.5), col = heat.colors(7)) #histogram of distances
abline(v=mean(distances(g, v="MuF",to=V(g))),col="blue", lty=3, lwd=2)

#Setting vertex color to something
clrs1=colorRampPalette(c("gray53","gray88","floralwhite","ghostwhite","white"))
cl1=clrs1(length(nodes$V1))
V(g)$color=cl1 #fill color
V(g)$frame.color="azure3" #frame color

#setting vertex color by degree
V(g)$color=ifelse(V(g)[degree(g)]>2,rainbow(length(degree(g)))[degree(g)],"gray")
V(g)$frame.color=V(g)$color #frame color

#scaling vertex size
V(g)$size=(degree(g)-min(degree(g)))/(max(degree(g))-min(degree(g)))*20+10 # scaled by degree

#scaling edge width
E(g)$width=graph_d$V3
E(g)$width=ifelse(log(E(g)$width)==0, .3,log(E(g)$width))

#coloring edges by width
ew=c(2.7,4.5)
E(g)$color=sapply(E(g)$width,
                  function(x) if(x>=ew[1] && x<=ew[2]) E(g)$color="cadetblue" else if(x>ew[2]) E(g)$color="maroon" else E(g)$color="gray")

#setting edge weight
E(g)$weight=E(g)$width

#removing low connectivity vertices
te=V(g)[degree(g)<2]
g1=delete_vertices(g,te)

#plotting
windows(48,32)
par(mai=c(0,0,0,0)) # set margins to fill plotting window
plot(g, vertex.label.dist=.2,asp=0,edge.curved=T)

tkplot(g,vertex.label.dist=0,asp=0,edge.curved=T)
l=tk_coords(tkp.id = 53)

#finding largest cliques
te=largest_cliques(g)
ye=unique(unlist(te))
ye=V(g)$name[ye]
g2=induced_subgraph(g,ye)

#finding cliques of particular size range
te=cliques(g,min = 7, max=8)
ye=unique(unlist(te)) # merging largest cliques
ye=V(g)$name[ye]
g2=induced_subgraph(g,ye) # getting clique subgraph

#setting vertex color by cliques

V(g)$color=ifelse(V(g)$name %in% ye,"orange","gray")
V(g)$frame.color=V(g)$color #frame color


tkplot(g,vertex.label.dist=.2,asp=0,edge.curved=T)

#community detection
com=cluster_edge_betweenness(g)
grp=com[4] # marking some community
grp=list(ye) # making group based on cliques found above

# png(width = 11, height = 8, units = 'in', res = 300)
par(mai=c(0,0,0,0))

plot(g, mark.groups = grp, vertex.size=V(g)$size/3,
     edge.curved=T, asp=0, vertex.label.cex=.3, layout=layout_with_lgl,
     edge.arrow.size=.3)
# dev.off()

#reordered graph
deg=degree(g)
nod=data.frame(V(g)$name[order(deg,decreasing = T)], sort(degree(g), decreasing = T))
colnames(nod)=c("node","degree")
nod=data.frame(nod$node,nod$degree)

g = graph_from_data_frame(graph_d, directed=TRUE, vertices = nod)

#replotting above as circular
n=15 # first circle
m=10 # second circle
t=seq(0,2*pi-2*pi/n,length=n)
x=rep(10*cos(t), m)
y=rep(10*sin(t), m)

t=seq(0,2*pi-2*pi/m,length=m)

z=40*cos(t)
z=as.vector(sapply(z,function(x) rep(x,n)))
x=x+z
z=40*sin(t)
z=as.vector(sapply(z,function(x) rep(x,n)))
y=y+z

V(g)$x=x
V(g)$y=y
V(g)$size=(degree(g)-min(degree(g)))/(max(degree(g))-min(degree(g)))*20+10 # scaled by degree
plot(g, vertex.label.dist=0, edge.curved=F, asp=1)
tkplot(g, vertex.label.dist=0, edge.curved=F, asp=1)

#can be done with unordered graph but best with reordered
#coloring nodes by rank of reaching in search algorithm bfs
te=bfs(g,"MuF",rank = T) #breadth first search
V(g)$color=ifelse(te$rank<100, rainbow(length(te$rank))[te$rank],"gray") # coloring 1st 100 nodes by rank
V(g)$frame.color=V(g)$color

#finding largest biconnected component
te=biconnected_components(g)
ye=max(unlist(lapply(te$components, function(x) length(x)))) # finding largest component
ye=grep(ye,lapply(te$components, function(x) length(x))) # finding component index with largest component
ye=unlist(te$components[ye]) # getting out that component
g1=induced_subgraph(g,ye) #create subgraph of largest component
tkplot(g1, vertex.label.dist=0, edge.curved=F, asp=1)

V(g)$color=ifelse(V(g)%in%ye, rainbow(length(V(g)))[V(g)],"grey") # color largest component in complete graph
V(g)$frame.color=V(g)$color

#sizing and coloring nodes by hub score or authority score
te=hub_score(g,scale = 1) #get hub score
te=authority_score(g,scale = 1) #get authority score instead

V(g)$size=te$vector*20+10
V(g)$color=ifelse(V(g)$size>15,rainbow(length(V(g)$size))[V(g)$size],"gray")
V(g)$frame.color=V(g)$color #frame color

#sizing and coloring nodes by betweenness
te=betweenness(g)
V(g)$size=(te-min(te))/(max(te)-min(te))*20+10
V(g)$color=ifelse(V(g)$size>10,rainbow(length(V(g)$size))[V(g)$size],"gray")
V(g)$frame.color=V(g)$color #frame color

tkplot(g, vertex.label.dist=0, edge.curved=T, asp=1, layout=layout_with_fr)
plot(g, vertex.label.dist=0, vertex.size=V(g)$size/3.5,
     edge.curved=T, asp=0, vertex.label.cex=.5, layout=coords,
     edge.arrow.width=.5)
