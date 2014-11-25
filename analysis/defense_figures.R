require(cowplot)
rm(list=ls())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols = c('gold', 'red')

draw_plot <- function(model, data, title)
{
  d <- cbind(data, predicted.w=predict(model, data))
  r.value <- cor(d$FEL.dN.dS, d$predicted.w, use="complete.obs")
  Rsq.label <- paste("R^2 == ", round(r.value^2, 2), sep='')
  print(Rsq.label)
  plot <- ggplot() +
  geom_point(data=d, aes(x=predicted.w, y=FEL.dN.dS, size=1/distance.to.224,     color=RSA.Multimer), alpha=0.75) +  
    scale_x_continuous(limits=c(0, 2.25), breaks=c(0, 1, 2)) +
    scale_colour_gradientn(colours = cols, name='    RSA') +
    scale_size(name='1 / distance\nto site 224', range = c(2, 6)) +
    ylim(0, 3.73) +
    xlab("Predicted dN/dS") +
    ylab("Observed dN/dS") +
    ggtitle(title) + 
    theme(plot.title = element_text(size=15, face="plain")) +
    geom_path(data=data.frame(x=c(0,2.25), y=c(0,2.25)), mapping=aes(x=x, y=y), color='black') +
    annotate("text", x=1.9, y=0.5, label=Rsq.label, parse=T, color='black')
  plot
}



d <- read.table('../manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)
d <- d[!is.na(d$pdb.4fnk), ]

m <- lm(FEL.dN.dS ~ Bush.99, data=d)
p1.1 <- draw_plot(m, d, "Bush '99")

m <- lm(FEL.dN.dS ~ RSA.Multimer + Bush.99, data=d)
p1.2 <- draw_plot(m, d, "Bush '99 + RSA")

m <- lm(FEL.dN.dS ~ Meyer.14, data=d)
p2.1 <- draw_plot(m, d, "Epitopes")

m <- lm(FEL.dN.dS ~ RSA.Multimer + Meyer.14, data=d)
p2.2 <- draw_plot(m, d, "Epitopes + RSA")

m <- lm(FEL.dN.dS ~ I(1 / distance.to.224), data=d)
p3.1 <- draw_plot(m, d, "1/Distance")

m <- lm(FEL.dN.dS ~ RSA.Multimer + I(1 / distance.to.224), data=d)
p3.2 <- draw_plot(m, d, "1/Distance + RSA")

p <- plot_grid(p1.1, p1.2, p2.1, p2.2, p3.1, p3.2, cols=2,
    labels = c("A", "B", "C", "D", "E", "F"))
    
ggsave("../analysis/combined_h3.pdf", p, width=11, height=12)
ggsave("../analysis/inverse_distance_h3.pdf", p3.1, width=5, height=4)
