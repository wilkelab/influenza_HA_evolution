rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)

b.edges <- read.table('epitope_data/b_cell_human_nonlinear.edges', head=T, sep='\t')
b.clusters <- epitopes$Pan.11[!is.na(epitopes$Protein)]
b.clusters[b.clusters=='-'] <- 'None'
b.vertices <- data.frame(site=1:length(b.clusters), ep=b.clusters)
b.vertices <- b.vertices[unique(c(b.edges$n1, b.edges$n2)), ]
b.vertices <- b.vertices[order(b.vertices$site), ]
b.network <- graph.data.frame(b.edges, vertices=b.vertices, directed=F)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('P' = mycols[1], 'None' = "#000000")

b <- simplify(b.network)
V(b)$color[V(b)$ep == 'P'] <- cbbPalette[1]
V(b)$color[V(b)$ep == 'None'] <- cbbPalette[2]
E(b)$weight <- seq(ecount(b))
E(b)$curved <- 0.1

pdf('analysis/b_human_nonlinear_network_11.pdf', height=7, width=7, useDingbats = F)
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
