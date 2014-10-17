rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)

#Nonhuman epitopes
b.linear.eps <- read.table('epitope_data/b_cell_nonhuman_linear.counts', head=T, stringsAsFactors = F)
b.nonlinear.eps <- read.table('epitope_data/b_cell_nonhuman_nonlinear.counts', head=T, stringsAsFactors = F)

binary.eps <- rep(NA, length(b.nonlinear.eps$count))
binary.eps[(b.nonlinear.eps$count + b.linear.eps$count) > 0] <- (b.nonlinear.eps$count + b.linear.eps$count)[(b.nonlinear.eps$count + b.linear.eps$count) > 0]
id <- epitopes$Pan.11[!is.na(epitopes$Protein)]
id[id=='-'] <- 'None'

df <- data.frame(index = epitopes$Protein[!is.na(epitopes$Protein)],
                 b.eps = b.linear.eps$count+b.nonlinear.eps$count, 
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('P' = mycols[1], 'A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'N' = mycols[6], 'M' = mycols[7], 'None' = "#000000")

p <- ggplot(aes(x=index, y=binary.eps, colour=id, fill=id), data=df) + geom_bar(stat='identity', width=0.1) + 
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + geom_tile(aes(x=index-0.01, y=rep(-0.4, length(index)), fill=id), height=0.4, width=0.75)
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550), expand=c(0,0))
p <- p + scale_y_continuous(breaks=seq(0, 20, 4), limits=c(-0.65, 18.1), expand=c(0,0))
p <- p + ylab('Count in Human B cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')

p <- p + theme(legend.position = c(0.8, 0.7),
               legend.title=element_blank())

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/b_nonhuman_epitopes_11.pdf', width=10, height=5)

#Human epitopes
b.linear.eps <- read.table('epitope_data/b_cell_human_linear.counts', head=T, stringsAsFactors = F)
b.nonlinear.eps <- read.table('epitope_data/b_cell_human_nonlinear.counts', head=T, stringsAsFactors = F)

binary.eps <- rep(NA, length(b.nonlinear.eps$count))
binary.eps[(b.nonlinear.eps$count + b.linear.eps$count) > 0] <- (b.nonlinear.eps$count + b.linear.eps$count)[(b.nonlinear.eps$count + b.linear.eps$count) > 0]

df <- data.frame(index = epitopes$Protein[!is.na(epitopes$Protein)],
                 b.eps = b.linear.eps$count+b.nonlinear.eps$count, 
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

p <- ggplot(aes(x=index, y=binary.eps, colour=id, fill=id), data=df) + geom_bar(stat='identity', width=0.1) + 
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + geom_tile(aes(x=index-0.01, y=rep(-0.4, length(index)), fill=id), height=0.4, width=0.75)
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550), expand=c(0,0))
p <- p + scale_y_continuous(breaks=seq(0, 20, 4), limits=c(-0.65, 18.1), expand=c(0,0))
p <- p + ylab('Count in Human B cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')

p <- p + theme(legend.position = c(0.8, 0.7),
               legend.title=element_blank())

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/b_human_epitopes_11.pdf', width=10, height=5)
