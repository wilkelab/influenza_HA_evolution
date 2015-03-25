rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)
set.seed(1)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)

b.edges <- read.table('epitope_data/b_cell_nonhuman_linear.edges', head=T, sep='\t')
b.clusters <- epitopes$Bush.99[!is.na(epitopes$Protein)]
b.clusters[b.clusters=='-'] <- 'None'
b.vertices <- data.frame(site=1:length(b.clusters), ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'N' = mycols[6], 'M' = mycols[7], 'None' = "#000000")

b <- simplify(b.network)
net <- graph.adjacency(get.adjacency(b.network),mode="undirected",weighted=TRUE,diag=FALSE)
V(net)$color[V(b)$ep == 'A'] <- cbbPalette[1]
V(net)$color[V(b)$ep == 'B'] <- cbbPalette[2]
V(net)$color[V(b)$ep == 'C'] <- cbbPalette[3]
V(net)$color[V(b)$ep == 'D'] <- cbbPalette[4]
V(net)$color[V(b)$ep == 'E'] <- cbbPalette[5]
V(net)$color[V(b)$ep == 'N'] <- cbbPalette[6]
V(net)$color[V(b)$ep == 'M'] <- cbbPalette[7]
V(net)$color[V(b)$ep == 'None'] <- cbbPalette[8]

l <- layout.fruchterman.reingold(net, niter=600, area=vcount(net)^2.4, repulserad=vcount(net)^2.9)

pdf('analysis/b_nonhuman_linear_network_99.pdf', height=7, width=7, useDingbats = F)
par(mar=c(0,0,0,0))
plot(net,
     layout=l,
     edge.width=0.1, 
     vertex.size = 3.5,
     vertex.frame.color= "white",
     vertex.label.color = "white",
     vertex.label.cex = 0.25,
     vertex.label.family = "sans",
     edge.color="black",
     edge.width=E(net)$weight
)

dev.off()
