rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
epitopes <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)
id <- rep('None', 550)
A <- c(34, 36, 50, 53, 54, 70, 292, 294, 386, 383, 305, 380, 398, 395, 382, 397, 365, 364, 363, 498, 307, 387, 393, 366, 394, 403, 384, 401, 405, 390, 334, 391, 404, 379)
B <- c(158,159,160,155,156,157,191,192,196,190,189,194,256,223,187,188)
C <- c(124,126,146,138,143,121,122,144,123,142,136,140,145,123,133,135,137,131)
D <- c(241,235,238,242,212,204,176,161,162,147,169,208,206,205,209,149,171,114,151,210,153,152,150,172,115,148,211,154,243,170,175,174,173,176)

id[A] <- '1'
id[B] <- '2'
id[C] <- '3'
id[D] <- '4'

require(latticeExtra)
mycols <- dput(ggplot2like(n = 7, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('1' = mycols[1], '2' = mycols[2], '3' = mycols[3], '4' = mycols[5], 'None' = "#000000")

#Human epitopes
b.nonlinear.eps <- read.table('epitope_data/b_cell_human_nonlinear.counts', head=T, stringsAsFactors = F)

binary.eps <- rep(NA, length(b.nonlinear.eps$count))
binary.eps[(b.nonlinear.eps$count) > 0] <- (b.nonlinear.eps$count)[(b.nonlinear.eps$count) > 0]

df <- data.frame(index = epitopes$Protein[!is.na(epitopes$Protein)],
                 b.eps = b.nonlinear.eps$count, 
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

p <- ggplot(aes(x=index, y=binary.eps, colour=id, fill=id), data=df) + geom_bar(stat='identity', width=0.1) + 
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + geom_tile(aes(x=index-0.01, y=rep(-0.2, length(index)), fill=id), height=0.3, width=0.75)
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550), expand=c(0,0))
p <- p + scale_y_continuous(breaks=seq(0, 20, 1), limits=c(-0.35, 5.1), expand=c(0,0))
p <- p + ylab('Count in Human B cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')

p <- p + theme(legend.position = c(0.8, 0.7),
               legend.title=element_blank())

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/b_human_epitopes_nonlinear_rename.pdf', width=10, height=5)

