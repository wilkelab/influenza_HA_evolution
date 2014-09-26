rm(list=ls())
seasonal.lines <- function(df) {
  require(ggplot2)
  require(grid)
  
  graphics.off()
  the_pointsize=18
  theme_set(theme_bw(base_size=the_pointsize))
  old_theme <- theme_update(panel.border=element_blank(),
                            axis.line=element_line(),
                            panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            axis.line=element_line())
  
  g <- ggplot(data=df, aes(x=year, y=y, color=id, fill=id)) + 
    geom_point(size=2.0) +
    geom_line(size=1.25) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.2, color=NA)
  
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Fraction of Sites in Epitopes')
  g <- g + xlab('Season')
  g <- g + scale_x_continuous(breaks=seq(1990, 2015, 5), limits=c(1992, 2015))
  g <- g + scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
  g <- g + theme(panel.border=element_blank(), axis.line=element_line())
  g <- g + theme(axis.title.x = element_text(size=24, vjust=-1))
  g <- g + theme(axis.text.x = element_text(size=24))
  g <- g + theme(axis.title.y = element_text(size=24, vjust=2))
  g <- g + theme(axis.text.y = element_text(size=24))
  g <- g + theme(axis.line = element_line(colour = 'black', size = 1))
  g <- g + theme(axis.ticks = element_line(colour = 'black', size = 1))
  g <- g + theme(plot.margin=unit(c(1.5, 1.5, 1.5, 1.5), "lines"))
  g <- g + theme(axis.ticks.margin = unit(0.25, "cm"))
  g <- g + theme(legend.position = "none")
  
  ggsave(g, file='~/Google Drive//Work/WilkeLab/flu_epitope_project/flu_chi_test/FracSiteInEps_P.pdf', width=10, height=10)
  return(g)
}

import.dat <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_chi_test/fraction_real_P.dat', sep='\t')
dat <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_chi_test/prob_data_P.dat', sep = ',')
dat <- t(dat)
colnames(dat) <- as.character(2:23)

alpha = 0.05
alpha = alpha/2
#me <- function(some.dat){qt(1 - alpha/2, length(some.dat) - 1) * sd(some.dat) / sqrt(length(some.dat))}
me <- function(some.dat){quantile(some.dat, c(alpha, 1-alpha))}
calc.ci <- apply(dat, 2, me)

means <- colMeans(dat)

colnames(import.dat) <- c('real_ep', 'rand_ep', 'ci')

plot.dat <- data.frame(id=c(rep('Real', length(calc.ci[1,])), rep('Random', length(calc.ci[1,]))), 
                       year=c(1993:2014, 1993:2014), 
                       y=c(import.dat$real_ep, means), 
                       ymax=c(import.dat$real_ep, calc.ci[2,]),
                       ymin=c(import.dat$real_ep, calc.ci[1,]))

seasonal.lines(plot.dat)