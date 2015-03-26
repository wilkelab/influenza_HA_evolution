rm(list=ls())
require(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = rev(gg_color_hue(3))

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/data_table/numbering_table_unix.csv', sep=',', head=T, stringsAsFactors = F)
d <- d[!is.na(d$pdb.4fnk), ]
r.value.1 <- summary(lm(SLAC.dN.All ~ RSA.Multimer, data=d))$r.squared
r.1 <- "RSA"

r.value.2 <- summary(lm(SLAC.dN.All ~ I(1/distance.to.224.all), data=d))$r.squared
r.2 <- "1 / Distance"

r.value.3 <- summary(lm(SLAC.dN.All ~ Bush.99, data=d))$r.squared
r.3 <- "Bush '99"

r.value.4 <- summary(lm(SLAC.dN.All ~ Meyer.14, data=d))$r.squared
r.4 <- "Epitopes"

r.value.5 <- summary(lm(SLAC.dN.All ~ I(1 / distance.to.224.all) + Bush.99, data=d))$r.squared
r.5 <- "1 / Distance + Bush '99"

r.value.6 <- summary(lm(SLAC.dN.All ~ I(1 / distance.to.224.all) + Meyer.14 , data=d))$r.squared
r.6 <- "1 / Distance + Epitopes"

r.value.7 <- summary(lm(SLAC.dN.All ~ RSA.Multimer + Bush.99, data=d))$r.squared
r.7 <- "RSA + Bush '99"

r.value.8 <- summary(lm(SLAC.dN.All ~ RSA.Multimer + Meyer.14 , data=d))$r.squared
r.8 <- "RSA + Epitopes"

r.value.9 <- summary(lm(SLAC.dN.All ~ RSA.Multimer + I(1 / distance.to.224.all), data=d))$r.squared
r.9 <- "RSA + 1 / Distance"

r.value.10 <- summary(lm(SLAC.dN.All ~ RSA.Multimer + I(1 / distance.to.224.all) + Bush.99, data=d))$r.squared
r.10 <- "RSA + 1 / Distance + Bush '99"

r.value.11 <- summary(lm(SLAC.dN.All ~ RSA.Multimer + I(1 / distance.to.224.all) + Meyer.14, data=d))$r.squared
r.11 <- "RSA + 1 / Distance + Epitopes"

rs <- c(r.value.1, r.value.4, r.value.2, r.value.3, r.value.8, r.value.7, r.value.6, r.value.5, r.value.9, r.value.11, r.value.10)
r.names <- c(r.1, r.4, r.2, r.3, r.8, r.7, r.6, r.5, r.9, r.11, r.10)

df <- data.frame(r.square = rs, names = factor(r.names, levels = r.names))

p <- ggplot(aes(x = names, y = r.square), data = df) +
  geom_bar(stat = 'identity', colour='darkgray', fill='darkgray') +
  ylab(expression(paste("Variance Explained (R"^"2", ')', sep=''))) +
  xlab("Predictor Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(0, 0.4))

ggsave("~/Google Drive/Data/influenza_HA_evolution/figures/r_squared_slac_all.pdf", p, width=7.5, height=7.5)