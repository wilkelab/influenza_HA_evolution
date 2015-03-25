rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)

df <- data.frame(index = epitopes$Protein[!is.na(epitopes$Protein)],
                 w = epitopes$FEL.dN.dS[!is.na(epitopes$Protein)]
)

p <- ggplot(aes(x=index, y=w), data=df) + geom_point()
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550))
p <- p + scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 4))
p <- p + ylab('dN/dS')
p <- p + xlab('Site in Mature Hemagglutinin')

p <- p + theme(legend.position = c(0.8, 0.7),
               legend.title=element_blank())

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/w.pdf', width=5, height=3.5)

