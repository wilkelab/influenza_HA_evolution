rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('A' = mycols[1], 'B' = mycols[2], 'C' = mycols[3], 'D' = mycols[4], 'E' = mycols[5], 'None' = "#000000")

setwd('~/Google Drive/Data/influenza_HA_evolution/')
ep.A<-c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168)
ep.B<-c(128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198)
ep.C<-c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312)
ep.D<-c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248)
ep.E<-c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)
ep.total<-c(ep.A, ep.B, ep.C, ep.D, ep.E)

get.data <- function(this.folder) {
  start <- this.folder
  dirs <- list.files(start)
  
  year <- c()
  rate <- c()
  index <- c()
  p <- c()
  ep <- c()
  
  for(i in dirs[1:length(dirs)]) {
    if(substr(i, 1, 5) == 'rates' & i != 'rates_combined.dat'){
      dat <- read.table(paste(start, i, sep=''), sep=',', header=T, stringsAsFactors=F)[17:566,]
      
      id <- rep('None', nrow(dat))
      id[ep.A] <- 'A'
      id[ep.B] <- 'B'
      id[ep.C] <- 'C'
      id[ep.D] <- 'D'
      id[ep.E] <- 'E'
      
      dat <- cbind(dat, index=1:nrow(dat), deparse.level=0)
      id <- id[dat$dN.dS > 1 & dat$p.value < 0.05]
      dat <- dat[dat$dN.dS > 1 & dat$p.value < 0.05, ]

      year <- append(year, rep(as.numeric(substr(i, 7, nchar(i) - 4)) + 1990, nrow(dat)))
      rate <- append(rate, dat$dN.dS)
      p <- append(p, dat$p.value)
      index <- append(index, dat$index)
      ep <- append(ep, id)
    }
  }
  
  return(data.frame(year=year, rate=rate, p=p, index=index, ep=ep, stringsAsFactors = FALSE))
}

df <- get.data('sequence_data/not_structure/omega/')
counts <- c()
years <- c()
id <- c()

for(i in unique(df$year)) {
  counts <- append(counts, as.vector(table(df$ep[df$year==i])))
  years <- append(years, rep(i, length(as.vector(table(df$ep[df$year==i])))))
  id <- append(id, unique(df$ep[df$year==i])[order(unique(df$ep[df$year==i]))])
}

plot.df <- data.frame(counts=counts, years=years, id=id)
plot.df <- plot.df[rev(order(plot.df$id)), ]

p <- ggplot(aes(x=years, y=counts, fill=factor(id)), data=plot.df) + geom_bar(stat='identity', position = "stack") +
  scale_colour_manual(values=cbbPalette) + 
  scale_fill_manual(values=cbbPalette)
p <- p + scale_x_continuous(breaks=seq(1991, 2014, 2), limits=c(1990, 2014), expand=c(0,0))
p <- p + scale_y_continuous(breaks=seq(0, 20, 4), limits=c(0, 18), expand=c(0,0))
p <- p + ylab('Count in T cell Epitopes')
p <- p + xlab('Site in Mature Hemagglutinin')

p <- p + theme(legend.position = c(0.95, 0.85),
               legend.title=element_blank())

ggsave(p, file='~/Google Drive/Data/influenza_HA_evolution/analysis/site_v_time.pdf', width=10, height=5)