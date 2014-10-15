rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
ep.A.skehel<-c(122, 144, 145, 137, 133, 143, 146, 132, 226)
ep.B.skehel<-c(155, 193, 188, 189, 186, 156, 157, 158, 126, 159, 160, 197, 198, 164)
ep.C.skehel<-c(18, 48, 50, 53, 54, 275, 91, 92, 275, 278)
ep.D.skehel<-c(173, 174, 172, 242, 244, 246, 248, 201, 205, 207, 208, 217, 220)
ep.E.skehel<-c(78, 81, 63, 83, 62)

ep.total<-c(ep.A, ep.B, ep.C, ep.D)
deem<-c(122, 124, 126, 131, 133, 137, 140, 142, 144, 145, 128, 155, 156, 157, 158, 159, 189, 192, 193, 197, 45, 50, 273, 275, 278, 312, 121, 172, 173, 201, 207, 219, 226, 227, 229, 246, 57, 62, 75, 78, 83, 92, 260, 262, 3, 5, 25, 33, 49, 106, 202, 222, 225, 271)

b.linear.eps <- (read.table('epitope_data/b_cell_nonhuman_linear.counts', head=T))
nonlinear.eps <- read.table('epitope_data/b_cell_nonhuman_nonlinear.counts', head=T)

binary.eps <- rep(NA, length(nonlinear.eps$count))
binary.eps[(nonlinear.eps$count + b.linear.eps$count) > 0] <- (nonlinear.eps$count + b.linear.eps$count)[(nonlinear.eps$count + b.linear.eps$count) > 0]
id <- rep('None', length(nonlinear.eps$count))
id[ep.A] <- 'A'
id[ep.B] <- 'B'
id[ep.C] <- 'C'
id[ep.D] <- 'D'
id[ep.E] <- 'E'

df <- data.frame(index = 1:550,
                 b.eps = b.linear.eps$count+nonlinear.eps$count, 
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

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

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/b_nonhuman_epitopes.pdf', width=10, height=5)

