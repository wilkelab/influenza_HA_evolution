require(cowplot)
rm(list=ls())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = rev(gg_color_hue(3))

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/data_table/numbering_table_unix.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ RSA.Multimer, data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p1.1 <- ggplot(data=d, aes(x=RSA.Multimer, y=FEL.dN.dS, size=1/distance.to.224), alpha=0.75) +
  geom_point() +
  stat_smooth(method='lm', se=F, color='red', size = 1) +
  scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  scale_size(name='1 / distance\nto site 224', range = c(2, 6)) +
  ylim(0, 3.73) +
  xlab("Relative Solvent Accessibility") +
  ylab("dN/dS") +
  annotate("text", x=0.15, y=3.5, label=paste("R^2 == ", round(r.value^2, 2), sep=''), parse=T, color='black')

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/data_table/numbering_table_unix.csv', sep=',', head=T, stringsAsFactors = F)

m <- lm(FEL.dN.dS ~ I(1/distance.to.224), data=d)
d <- cbind(d, predicted.w=predict(m, d))
r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
p1.2 <- ggplot(data=d, aes(x=(1/distance.to.224), y=FEL.dN.dS, color=RSA.Multimer), alpha=0.75) +
  geom_point() +
  stat_smooth(method='lm', se=F, color='black', size = 1) +
  scale_x_continuous(limits=c(0, 0.3), breaks=seq(0, 0.25, by=0.05)) +
  scale_colour_gradientn(colours = c('gold', 'red'), name='RSA') +
  ylim(0, 3.73) +
  xlab("Inverse Distance to Site 224 (1/Angstrom)") +
  ylab("dN/dS") +
  annotate("text", x=0.04, y=3.5, label=paste("R^2 == ", round(r.value^2, 2), sep=''), parse=T, color='black')

p <- plot_grid(p1.1, p1.2, cols=2)
p <- p + draw_plot_label(c("A", "B"), c(0, 1/2), c(1, 1))
ggsave("~/Google Drive/Data/influenza_HA_evolution/analysis/rsa_inversedistance.pdf", p, width=12, height=5)