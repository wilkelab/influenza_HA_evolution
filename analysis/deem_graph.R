rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
deem<-c(122, 124, 126, 131, 133, 137, 140, 142, 144, 145, 128, 155, 156, 157, 158, 159, 189, 192, 193, 197, 45, 50, 273, 275, 278, 312, 121, 172, 173, 201, 207, 219, 226, 227, 229, 246, 57, 62, 75, 78, 83, 92, 260, 262, 3, 5, 25, 33, 49, 106, 202, 222, 225, 271)

b.edges <- read.table('epitope_data/b_edges.dat', head=T, sep=',')
b.clusters <- rep('None', 549)
b.clusters[deem] <- 'Deem'
b.vertices <- data.frame(site=1:549, ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('Deem' = mycols[1], 'None' = "#000000")

b <- simplify(b.network)
V(b)$color[V(b)$ep == 'Deem'] <- cbbPalette[1]
V(b)$color[V(b)$ep == 'None'] <- cbbPalette[2]
E(b)$weight <- seq(ecount(b))
E(b)$curved <- 0.1

pdf('analysis/deem_network.pdf', height=7, width=7, useDingbats = F)
par(mar=c(0,0,0,0))
plot(b,
     #layout=layout.circle,
     edge.width=0.2, 
     vertex.size = 3,
     vertex.frame.color= "white",
     vertex.label.color = "white",
     vertex.label.cex = 0.2,
     vertex.label.family = "sans",
     edge.color="black")

dev.off()

