rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)
set.seed(3)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
b.edges <- read.table('sequence_data/structure/4fnk_monomer.edges', head=T, sep=',')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)
b.clusters <- epitopes$Bush.99[!is.na(epitopes$Protein)]
b.clusters[b.clusters=='-'] <- 'None'
b.vertices <- data.frame(site=1:length(b.clusters), ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

b <- simplify(b.network)
net <- graph.adjacency(get.adjacency(b.network),mode="undirected",weighted=TRUE,diag=FALSE)
V(net)$color[V(b)$ep == 'A'] <- cbbPalette[1]
V(net)$color[V(b)$ep == 'B'] <- cbbPalette[2]
V(net)$color[V(b)$ep == 'C'] <- cbbPalette[3]
V(net)$color[V(b)$ep == 'D'] <- cbbPalette[4]
V(net)$color[V(b)$ep == 'E'] <- cbbPalette[5]
V(net)$color[V(b)$ep == 'None'] <- cbbPalette[6]

l <- layout.fruchterman.reingold(net, niter=600, area=vcount(net)^2.3, repulserad=vcount(net)^3.0)

pdf('analysis/b_nonhuman_nonlinear_distance_network_99.pdf', height=7, width=7, useDingbats = F)
par(mar=c(0,0,0,0))
plot(net,
     layout=l,
     edge.width=0.1, 
     vertex.size = 6,
     vertex.frame.color= "white",
     vertex.label.color = "white",
     vertex.label.cex = 0.45,
     vertex.label.family = "sans",
     edge.color="black",
     edge.width=E(net)$weight
)

dev.off()
