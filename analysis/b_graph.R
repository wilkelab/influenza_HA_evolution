rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
ep.A<-c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168)
ep.B<-c(128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198)
ep.C<-c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312)
ep.D<-c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248)
ep.E<-c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)
ep.total<-c(ep.A, ep.B, ep.C, ep.D, ep.E)

b.edges <- read.table('epitope_data/b_cell_nonhuman_linear.edges', head=T, sep='\t')
b.clusters <- rep('None', 550)
b.clusters[ep.A] <- 'A'
b.clusters[ep.B] <- 'B'
b.clusters[ep.C] <- 'C'
b.clusters[ep.D] <- 'D'
b.clusters[ep.E] <- 'E'
b.vertices <- data.frame(site=1:550, ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

b <- simplify(b.network)
V(b)$color[V(b)$ep == 'A'] <- cbbPalette[1]
V(b)$color[V(b)$ep == 'B'] <- cbbPalette[2]
V(b)$color[V(b)$ep == 'C'] <- cbbPalette[3]
V(b)$color[V(b)$ep == 'D'] <- cbbPalette[4]
V(b)$color[V(b)$ep == 'E'] <- cbbPalette[5]
V(b)$color[V(b)$ep == 'None'] <- cbbPalette[6]
E(b)$weight <- seq(ecount(b))
E(b)$curved <- 0.1

pdf('analysis/b_nonhuman_linear_network.pdf', height=7, width=7, useDingbats = F)
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
