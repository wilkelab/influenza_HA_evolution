require(cowplot)
rm(list=ls())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = rev(gg_color_hue(3))

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ RSA.Multimer, data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p1.1 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1/distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') +
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ RSA.Multimer + Bush.99, data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p1.2 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1/distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') +
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ I(1/distance.to.224), data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p2.1 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1 / distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') + 
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ I(1 / distance.to.224) + Bush.99, data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p2.2 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1 / distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') + 
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)
#d$Meyer.14[d$Meyer.14 == "N"] <- "-"

m <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224) + Meyer.14 , data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p3.1 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1/distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') +
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224) + Bush.99, data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p3.2 <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, color=1/distance.to.224, size=RSA.Multimer), alpha=0.75) +
  scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
  scale_colour_gradientn(colours = cols, name='1 / distance to\nresidue 224') +
  scale_size(name='    RSA', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Predicted dN/dS") +
  ylab("Observed dN/dS") +
  geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
  annotate("text", x=1.9, y=0.5, label=paste("r = ", round(r.value, 3), "***", sep=""), color='black')

p <- plot_grid(p1.1, p1.2, p2.1, p2.2, p3.1, p3.2, cols=2)
p <- p + draw_plot_label(c("A", "B", "C", "D", "E", "F"), c(0, 1/2, 0, 1/2, 0, 1/2), c(1, 1, 2/3, 2/3, 1/3, 1/3))
ggsave("~/Google Drive/Data/influenza_HA_evolution/analysis/combined_h3.pdf", p, width=11, height=12)