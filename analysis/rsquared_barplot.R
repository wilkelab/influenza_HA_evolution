rm(list=ls())
require(cowplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = rev(gg_color_hue(3))

d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

m1 <- lm(FEL.dN.dS ~ RSA.Multimer, data=d)
d1 <- cbind(d, predicted.w=predict(m1, d))
r.value.1 <- cor(d1$FEL.dN.dS, d1$predicted.w, use="complete.obs")
r.1 <- "RSA"

m2 <- lm(FEL.dN.dS ~ I(1/distance.to.224), data=d)
d2 <- cbind(d, predicted.w=predict(m2, d))
r.value.2 <- cor(d2$FEL.dN.dS, d2$predicted.w, use="complete.obs")
r.2 <- "1 / Distance"

m3 <- lm(FEL.dN.dS ~ Bush.99, data=d)
d3 <- cbind(d, predicted.w=predict(m3, d))
r.value.3 <- cor(d3$FEL.dN.dS, d3$predicted.w, use="complete.obs")
r.3 <- "Bush '99"

m4 <- lm(FEL.dN.dS ~ Meyer.14, data=d)
d4 <- cbind(d, predicted.w=predict(m4, d))
r.value.4 <- cor(d4$FEL.dN.dS, d4$predicted.w, use="complete.obs")
r.4 <- "Epitopes"

m5 <- lm(FEL.dN.dS ~ I(1 / distance.to.224) + Bush.99, data=d)
d5 <- cbind(d, predicted.w=predict(m5, d))
r.value.5 <- cor(d5$FEL.dN.dS, d5$predicted.w, use="complete.obs")
r.5 <- "1 / Distance + Bush '99"

m6 <- lm(FEL.dN.dS ~ I(1 / distance.to.224) + Meyer.14 , data=d)
d6 <- cbind(d, predicted.w=predict(m6, d))
r.value.6 <- cor(d6$FEL.dN.dS, d6$predicted.w, use="complete.obs")
r.6 <- "1 / Distance + Epitopes"

m7 <- lm(FEL.dN.dS ~ RSA.Multimer + Bush.99, data=d)
d7 <- cbind(d, predicted.w=predict(m7, d))
r.value.7 <- cor(d7$FEL.dN.dS, d7$predicted.w, use="complete.obs")
r.7 <- "RSA + Bush '99"

m8 <- lm(FEL.dN.dS ~ RSA.Multimer + Meyer.14 , data=d)
d8 <- cbind(d, predicted.w=predict(m8, d))
r.value.8 <- cor(d8$FEL.dN.dS, d8$predicted.w, use="complete.obs")
r.8 <- "RSA + Epitopes"

m9 <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224), data=d)
d9 <- cbind(d, predicted.w=predict(m9, d))
r.value.9 <- cor(d9$FEL.dN.dS, d9$predicted.w, use="complete.obs")
r.9 <- "RSA + 1 / Distance"

m10 <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224) + Bush.99, data=d)
d10 <- cbind(d, predicted.w=predict(m10, d))
r.value.10 <- cor(d10$FEL.dN.dS, d10$predicted.w, use="complete.obs")
r.10 <- "RSA + 1 / Distance + Bush '99"

m11 <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224) + Meyer.14, data=d)
d11 <- cbind(d, predicted.w=predict(m11, d))
r.value.11 <- cor(d11$FEL.dN.dS, d11$predicted.w, use="complete.obs")
r.11 <- "RSA + 1 / Distance + Epitopes"

rs <- c(r.value.1, r.value.2, r.value.3, r.value.4, r.value.8, r.value.7, r.value.6, r.value.5, r.value.9, r.value.11, r.value.10)^2
r.names <- c(r.1, r.2, r.3, r.4, r.8, r.7, r.6, r.5, r.9, r.11, r.10)

df <- data.frame(r.square = rs, names = factor(r.names, levels = r.names))

p <- ggplot(aes(x = names, y = r.square), data = df) +
  geom_bar(stat = 'identity', colour='darkgray', fill='darkgray') +
  ylab(expression(paste("Variance Explained (R"^"2", ')', sep=''))) +
  xlab("Predictor Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(limits = c(0, 0.4))

ggsave("~/Google Drive/Data/influenza_HA_evolution/analysis/r_squared.pdf", p, width=7.5, height=7.5)