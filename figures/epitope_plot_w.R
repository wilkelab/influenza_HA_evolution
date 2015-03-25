rm(list = ls())

ep.A<-c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168)
ep.B<-c(128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198)
ep.C<-c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312)
ep.D<-c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248)
ep.E<-c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)
ep.total<-c(ep.A, ep.B, ep.C, ep.D, ep.E)
ep.exp<-c(18,  20,  34,  37,  38,  54,  66,  67,  98, 105, 106, 107, 108, 110, 115, 117, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130, 131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 155, 156, 157, 158, 159, 160, 173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 222, 225, 226, 259, 260, 262, 276, 278, 289, 291, 318, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 347, 348, 349, 350, 363, 364, 366, 367, 368, 370, 371, 374, 375, 377, 378, 379, 381, 382, 385, 387, 388, 389, 482)
deem<-c(122, 124, 126, 131, 133, 137, 140, 142, 144, 145, 128, 155, 156, 157, 158, 159, 189, 192, 193, 197, 45, 50, 273, 275, 278, 312, 121, 172, 173, 201, 207, 219, 226, 227, 229, 246, 57, 62, 75, 78, 83, 92, 260, 262, 3, 5, 25, 33, 49, 106, 202, 222, 225, 271)

require(latticeExtra)
mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('Background' = "#000000", 'A' = mycols[1], 'B' = mycols[2], 'D' = mycols[3], 'C' = mycols[4], 'E' = mycols[5], 'Chi' = mycols[1], 'Deem' = mycols[1], 'B-Cell' = mycols[1])

w.lines <- function(df, name) {
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
  
  g <- ggplot(df, aes(x=num, y=w, fill=id, color=id)) + 
    geom_bar(stat='identity', width=1) + 
    geom_point() + 
    scale_colour_manual(values=cbbPalette) + 
    scale_fill_manual(values=cbbPalette)
  g <- g + geom_hline(yintercept=1)
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('dN/dS')
  g <- g + xlab('Index')
  g <- g + scale_x_continuous(breaks=seq(0, 1000, 100), limits=c(-20, 565))
  #  g <- g + scale_y_continuous(breaks=seq(500, 1000, 100), limits=c(550, 900))
  g <- g + theme(panel.border=element_blank(), axis.line=element_line())
  g <- g + theme(axis.title.x = element_text(size=24, vjust=-1))
  g <- g + theme(axis.text.x = element_text(size=24))
  g <- g + theme(axis.title.y = element_text(size=24, vjust=2))
  g <- g + theme(axis.text.y = element_text(size=24))
  g <- g + theme(axis.line = element_line(colour = 'black', size = 1))
  g <- g + theme(axis.ticks = element_line(colour = 'black', size = 1))
  g <- g + theme(plot.margin=unit(c(1.5, 1.5, 1.5, 1.5), "lines"))
  g <- g + theme(axis.ticks.margin = unit(0.25, "cm"))
  #g <- g + theme(legend.position = "none")
  
  ggsave(g, file=paste('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/', name, sep=''), width=12, height=10)
  return(g)
}

prob.w.91 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/1/result_after.dat', skip=1, sep=',')
prob.w.92 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/2/result_after.dat', skip=1, sep=',')
prob.w.93 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/3/result_after.dat', skip=1, sep=',')
prob.w.94 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/4/result_after.dat', skip=1, sep=',')
prob.w.95 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/5/result_after.dat', skip=1, sep=',')
prob.w.96 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/6/result_after.dat', skip=1, sep=',')
prob.w.97 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/7/result_after.dat', skip=1, sep=',')
prob.w.98 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/8/result_after.dat', skip=1, sep=',')
prob.w.99 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/9/result_after.dat', skip=1, sep=',')
prob.w.00 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/10/result_after.dat', skip=1, sep=',')
prob.w.01 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/11/result_after.dat', skip=1, sep=',')
prob.w.02 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/12/result_after.dat', skip=1, sep=',')
prob.w.03 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/13/result_after.dat', skip=1, sep=',')
prob.w.04 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/14/result_after.dat', skip=1, sep=',')
#prob.w.05 <- read.table('15/result_after.dat', skip=1, sep=',')
prob.w.06 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/16/result_after.dat', skip=1, sep=',')
prob.w.07 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/17/result_after.dat', skip=1, sep=',')
prob.w.08 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/18/result_after.dat', skip=1, sep=',')
prob.w.09 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/19/result_after.dat', skip=1, sep=',')
prob.w.10 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/20/result_after.dat', skip=1, sep=',')
prob.w.11 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/21/result_after.dat', skip=1, sep=',')
prob.w.12 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/22/result_after.dat', skip=1, sep=',')
prob.w.13 <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/23/result_after.dat', skip=1, sep=',')

source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/1/LnLik.dat')
w.91 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/2/LnLik.dat')
w.92 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/3/LnLik.dat')
w.93 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/4/LnLik.dat')
w.94 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/5/LnLik.dat')
w.95 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/6/LnLik.dat')
w.96 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/7/LnLik.dat')
w.97 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/8/LnLik.dat')
w.98 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/9/LnLik.dat')
w.99 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/10/LnLik.dat')
w.00 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/11/LnLik.dat')
w.01 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/12/LnLik.dat')
w.02 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/13/LnLik.dat')
w.03 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/14/LnLik.dat')
w.04 <- c(w_1, w_2, w_3, w_4)
#source('15/LnLik.dat')
#w.05 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/16/LnLik.dat')
w.06 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/17/LnLik.dat')
w.07 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/18/LnLik.dat')
w.08 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/19/LnLik.dat')
w.09 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/20/LnLik.dat')
w.10 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/21/LnLik.dat')
w.11 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/22/LnLik.dat')
w.12 <- c(w_1, w_2, w_3, w_4)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/23/LnLik.dat')
w.13 <- c(w_1, w_2, w_3, w_4)

total.w.91 <- colSums(w.91*prob.w.91)
total.w.92 <- colSums(w.92*prob.w.92)
total.w.93 <- colSums(w.93*prob.w.93)
total.w.94 <- colSums(w.94*prob.w.94)
total.w.95 <- colSums(w.95*prob.w.95)
total.w.96 <- colSums(w.96*prob.w.96)
total.w.97 <- colSums(w.97*prob.w.97)
total.w.98 <- colSums(w.98*prob.w.98)
total.w.99 <- colSums(w.99*prob.w.99)
total.w.00 <- colSums(w.00*prob.w.00)
total.w.01 <- colSums(w.01*prob.w.01)
total.w.02 <- colSums(w.02*prob.w.02)
total.w.03 <- colSums(w.03*prob.w.03)
total.w.04 <- colSums(w.04*prob.w.04)
#total.w.05 <- colSums(w.05*prob.w.05)
total.w.06 <- colSums(w.06*prob.w.06)
total.w.07 <- colSums(w.07*prob.w.07)
total.w.08 <- colSums(w.08*prob.w.08)
total.w.09 <- colSums(w.09*prob.w.09)
total.w.10 <- colSums(w.10*prob.w.10)
total.w.11 <- colSums(w.11*prob.w.11)
total.w.12 <- colSums(w.12*prob.w.12)
total.w.13 <- colSums(w.13*prob.w.13)

total.w <- data.frame(id=c(rep('91', 566), rep('92', 566), rep('93', 566), 
                           rep('94', 566), rep('95', 566), rep('96', 566), 
                           rep('97', 566), rep('98', 566), rep('99', 566),
                           rep('00', 566), rep('01', 566), rep('02', 566), 
                           rep('03', 566), rep('04', 566), #rep('05', 565), 
                           rep('06', 566), rep('07', 566), rep('08', 566),
                           rep('09', 566), rep('10', 566), rep('11', 566),
                           rep('12', 566), rep('13', 566)),
                      w=c(total.w.91, total.w.92, total.w.93, total.w.94, 
                          total.w.95, total.w.96, total.w.97, total.w.98,
                          total.w.99, total.w.00, total.w.01, total.w.02,
                          total.w.03, total.w.04, #total.w.05, 
                          total.w.06,
                          total.w.07, total.w.08, total.w.09, total.w.10,
                          total.w.11, total.w.12, total.w.13)
                      )
                      
prob.w.combined <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/combined/result_after.dat', sep=',', skip=1)
source('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_w_test/combined/LnLik.dat')
w.combined <- c(w_1,w_2,w_3,w_4)
#Need to fix the site number for the mature protein... Site Q17 in the transcript is amino acid 1 in the mature
total.w.combined <- data.frame(site=-15:(length(prob.w.combined)-16), w=colSums(w.combined*prob.w.combined))

tmp.data <- read.table('~/Google Drive/Work/WilkeLab/flu_epitope_project/flu_chi_test/positive_sites_all.dat')
site.counts <- table(tmp.data[,1])
site.counts <- as.numeric(names(site.counts[site.counts>1]))
positive.sites <- sort(unique(tmp.data)[, 1])

id = rep('Background', length(total.w.combined$site))
#id[total.w.combined$site %in% ep.A] = 'A'
#id[total.w.combined$site %in% ep.B] = 'B'
#id[total.w.combined$site %in% ep.C] = 'C'
#id[total.w.combined$site %in% ep.D] = 'D'
#id[total.w.combined$site %in% ep.E] = 'E'
#id[total.w.combined$site %in% site.counts] = 'Chi'
#id[total.w.combined$site %in% deem] = 'Deem'
id[total.w.combined$site %in% ep.exp] = 'B-Cell'

ep.w <- total.w.combined[total.w.combined$site %in% ep.total,]$w
ep.exp.w <- total.w.combined[total.w.combined$site %in% ep.exp,]$w
positive.w <- total.w.combined[total.w.combined$site %in% positive.sites,]$w
multiple.w <- total.w.combined[total.w.combined$site %in% site.counts,]$w
deem.w <- total.w.combined[total.w.combined$site %in% deem,]$w

plot.data <- data.frame(num=total.w.combined$site, w=total.w.combined$w, id=id)

w.lines(plot.data, name='bcell_plot.pdf')                      