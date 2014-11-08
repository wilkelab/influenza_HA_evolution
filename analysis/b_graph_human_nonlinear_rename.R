rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
A <- c(34, 36, 50, 53, 54, 70, 292, 294, 386, 383, 305, 380, 398, 395, 382, 397, 365, 364, 363, 498, 307, 387, 393, 366, 394, 403, 384, 401, 405, 390, 334, 391, 404, 379)
B <- c(158,159,160,155,156,157,191,192,196,190,189,194,256,223,187,188)
C <- c(124,126,146,138,143,121,122,144,123,142,136,140,145,123,133,135,137,131)
D <- c(241,235,238,242,212,204,176,161,162,147,169,208,206,205,209,149,171,114,151,210,153,152,150,172,115,148,211,154,243,170,175,174,173,176)
b.edges <- read.table('epitope_data/b_cell_human_nonlinear.edges', head=T, sep='\t')
b.clusters <- rep('None', 550)
b.clusters[A] <- '1'
b.clusters[B] <- '2'
b.clusters[C] <- '3'
b.clusters[D] <- '4'
b.vertices <- data.frame(site=1:length(b.clusters), ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('1' = mycols[1], '2' = mycols[2], '3' = mycols[3], '4' = mycols[5], 'None' = "#000000")

b <- simplify(b.network)
net <- graph.adjacency(get.adjacency(b.network),mode="undirected",weighted=TRUE,diag=FALSE)
V(net)$color[V(b)$ep == '1'] <- cbbPalette[1]
V(net)$color[V(b)$ep == '2'] <- cbbPalette[2]
V(net)$color[V(b)$ep == '3'] <- cbbPalette[3]
V(net)$color[V(b)$ep == '4'] <- cbbPalette[4]
V(net)$color[V(b)$ep == 'None'] <- cbbPalette[5]

l <- layout.fruchterman.reingold(net, niter=600, area=vcount(net)^2.3, repulserad=vcount(net)^3.0)

pdf('analysis/b_human_nonlinear_network_rename.pdf', height=7, width=7, useDingbats = F)
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
