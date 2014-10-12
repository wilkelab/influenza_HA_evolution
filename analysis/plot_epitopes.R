rm(list = ls())

ep.A<-c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168)
ep.B<-c(128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198)
ep.C<-c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312)
ep.D<-c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248)
ep.E<-c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)
ep.total<-c(ep.A, ep.B, ep.C, ep.D, ep.E)
ep.exp<-c(18,  20,  34,  37,  38,  54,  66,  67,  98, 105, 106, 107, 108, 110, 115, 117, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130, 131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 155, 156, 157, 158, 159, 160, 173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 222, 225, 226, 259, 260, 262, 276, 278, 289, 291, 318, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 347, 348, 349, 350, 363, 364, 366, 367, 368, 370, 371, 374, 375, 377, 378, 379, 381, 382, 385, 387, 388, 389, 482)
deem<-c(122, 124, 126, 131, 133, 137, 140, 142, 144, 145, 128, 155, 156, 157, 158, 159, 189, 192, 193, 197, 45, 50, 273, 275, 278, 312, 121, 172, 173, 201, 207, 219, 226, 227, 229, 246, 57, 62, 75, 78, 83, 92, 260, 262, 3, 5, 25, 33, 49, 106, 202, 222, 225, 271)

t.linear.eps <- (read.table('../epitope_data/t_linear.counts', head=T))[17:566, ]
b.linear.eps <- (read.table('../epitope_data/b_linear.counts', head=T))[17:566, ]
nonlinear.eps <- read.table('../epitope_data/nonlinear_eps.counts', head=T)
map <- t(read.table('../sequence_data/not_structure/combined/sequences/map.txt', sep=','))

#The epitopes already have the first 16 amino acids cut off for numbering
#linear.eps.short <- linear.eps$ep_counts[!is.na(as.vector(map[17:566,1]))]
#nonlinear.eps.short <- nonlinear.eps$ep_counts[!is.na(as.vector(map[17:566,1]))]

binary.eps <- rep(0, length(nonlinear.eps$ep_counts))
binary.eps[(nonlinear.eps$ep_counts + b.linear.eps) > 0 ] <- 1
id <- rep('None', length(nonlinear.eps$ep_counts))
id[ep.A] <- 'A'
id[ep.B] <- 'B'
id[ep.C] <- 'C'
id[ep.D] <- 'D'
id[ep.E] <- 'E'
#id[deem] <- 'Deem'

df <- data.frame(index = 1:550,
                 b.eps = b.linear.eps+nonlinear.eps$ep_counts, 
                 t.eps = t.linear.eps,
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'D' = mycols[3], 'C' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

p <- ggplot(aes(x=index, y=b.eps, colour=id, fill=id), data=df) + geom_bar(stat='identity', width=1) + geom_point(size=3) + 
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550))
p <- p + scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 5))
p <- p + ylab('Count in B cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')
p <- p + theme_bw()
p <- p + theme(panel.border=element_blank())
p <- p + theme(panel.border=element_blank(), 
               axis.line=element_line(),
               axis.title.x = element_text(size=32, vjust=0),
               axis.text.x = element_text(size=24),
               axis.title.y = element_text(size=32, vjust=1.5),
               axis.text.y = element_text(size=24),
               axis.line = element_line(colour = 'black', size = 1),
               axis.ticks = element_line(colour = 'black', size = 1),
               plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
               axis.ticks.margin = unit(0.1, "cm")
)
p <- p + theme(legend.position = c(0.8, 0.8),
               legend.title=element_blank(), 
               legend.key = element_blank(), 
               legend.text=element_text(size=24),
               legend.key.size = unit(1, "cm"))

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/b_epitopes.pdf', width=20, height=10, useDingbats=FALSE)

##T cell epitopes
binary.eps <- rep(0, length(nonlinear.eps$ep_counts))
binary.eps[(t.linear.eps) > 0 ] <- 1
id <- rep('None', length(t.linear.eps))
id[ep.A] <- 'A'
id[ep.B] <- 'B'
id[ep.C] <- 'C'
id[ep.D] <- 'D'
id[ep.E] <- 'E'
#id[deem] <- 'Deem'

df <- data.frame(index = 1:550,
                 b.eps = b.linear.eps+nonlinear.eps$ep_counts, 
                 t.eps = t.linear.eps,
                 binary.eps = binary.eps,
                 id=id
)

df <- df[rev(order(df$id)), ]

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'D' = mycols[3], 'C' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

p <- ggplot(aes(x=index, y=t.eps, colour=id, fill=id), data=df) + geom_bar(stat='identity', width=1) + geom_point(size=3) + 
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + scale_x_continuous(breaks=seq(0, 550, 100), limits=c(0, 550))
#p <- p + scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 5))
p <- p + ylab('Count in T cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')
p <- p + theme_bw()
p <- p + theme(panel.border=element_blank())
p <- p + theme(panel.border=element_blank(), 
               axis.line=element_line(),
               axis.title.x = element_text(size=32, vjust=0),
               axis.text.x = element_text(size=24),
               axis.title.y = element_text(size=32, vjust=1.5),
               axis.text.y = element_text(size=24),
               axis.line = element_line(colour = 'black', size = 1),
               axis.ticks = element_line(colour = 'black', size = 1),
               plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
               axis.ticks.margin = unit(0.1, "cm")
)
p <- p + theme(legend.position = c(0.8, 0.8),
               legend.title=element_blank(), 
               legend.key = element_blank(), 
               legend.text=element_text(size=24),
               legend.key.size = unit(1, "cm"))

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/t_epitopes.pdf', width=20, height=10, useDingbats=FALSE)