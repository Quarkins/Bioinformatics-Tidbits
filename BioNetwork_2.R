#A quick test to build a network for bioinformatics server

library("igraph")
library("network")
library("sna")
library("ndtv")
library("RColorBrewer")

#read in data
data = read.csv(file="~/Downloads/Data Science collaborations - Sheet1.csv",header=TRUE,stringsAsFactors = FALSE)

node_names = unique(data[,"Collaborator.group"])
node_ids = seq(1,length(node_names),by=1)
theme_names = unique(data[,c("Collaborator.theme","Collaborator.group")])

node_themes = vector()
#Match node_names to theme
for(i in 1:length(node_names)){
  node_themes[i] = theme_names[i,1]
}

#Add node for bioinformatics
node_names = c(node_names,"Bioinformatics")
node_ids = c(node_ids,length(node_ids)+1)
node_themes = c(node_themes,"Bioinformatics")

#Create node dataframe
nodes = data.frame(id=node_ids,name=node_names,theme=node_themes)

#Create link data frame
links <- data.frame(from=integer(),to=integer(),weight=integer())
for(i in 1:nrow(data)){
  from =  which(data[i,1] == node_names)
  to = which(data[i,4] == node_names)
  weight = 1
  links[i,] = c(from,to,weight)
}

#Force only one link per node
links2 = links[!duplicated(links[,1:2]),]

#Calculate the new weights (based on the number of connections)
node_size = count(links,"to")
nodes_si = c((node_size$freq)^2+20,21) #Add the one for bioninformatics

net <- graph.data.frame(links2,nodes,directed=FALSE)

#Color nodes by group

#pal3 <- colorRampPalette(brewer.pal(12,"Spectral"))(length(node_ids))
pal3 <- brewer.pal(length(unique(node_themes)),"Set1")
V(net)$color <- pal3[as.integer(factor(node_themes))]
mcri_pal = c("#00ADEF","#8DC63F","#00B7C6","#F47920","#7A52C7","#EC008C")
V(net)$color <- mcri_pal[as.integer(factor(node_themes))]

#Get degree
#deg = igraph::degree(net,mode="all")

V(net)$size <- nodes_si
E(net)$width <- E(net)$weight #Edge width defined by connection

#Layout
layoutVertexByAttr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph))
  V(g)$name <- as.character(1:vcount(g))
  E(g)$weight <- 1
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  l <- layout(g,weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}

l2 = layout.fruchterman.reingold(net)
l3 = layoutVertexByAttr(net, V(net)$color, cluster.strength=100,
                        layout=layout.fruchterman.reingold)
pdf("connections2.pdf",width=12,height=12)
plot(net,vertex.label=NA,layout=l3,
     edge.color=adjustcolor("#00ADEF", alpha.f = .5),
     vertex.color=adjustcolor(V(net)$color,alpha.f=.9))
legend("bottomright",legend=levels(factor(node_themes)),
       pch=21,col="#777777",
       pt.bg=mcri_pal, pt.cex=1, cex=1, bty="n", ncol=1,inset=.05)
dev.off()
