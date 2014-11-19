rm(list = ls())
d <- read.table('~/Google Drive/Data/influenza_HA_evolution/manuscript/numbering_table.csv', sep=',', head=T, stringsAsFactors = F)

dist.non_epitope.99 <- d$distance.to.224[d$Bush.99 == '-']
dist.epitope.99 <- d$distance.to.224[d$Bush.99 != '-']
dist.non_epitope.rename <- d$distance.to.224[d$Meyer.14 == '-']
dist.epitope.rename <- d$distance.to.224[d$Meyer.14 != '-']

rsa.non_epitope.99 <- d$RSA.Multimer[d$Bush.99 == '-']
rsa.epitope.99 <- d$RSA.Multimer[d$Bush.99 != '-']
rsa.non_epitope.rename <- d$RSA.Multimer[d$Meyer.14 == '-']
rsa.epitope.rename <- d$RSA.Multimer[d$Meyer.14 != '-']

df <- data.frame(distances = c(dist.epitope.99, dist.non_epitope.99, dist.epitope.rename, dist.non_epitope.rename),
                 rsa = c(rsa.epitope.99, rsa.non_epitope.99, rsa.epitope.rename, rsa.non_epitope.rename), 
                 id = c(rep("Epitopes Bush 1999", length(dist.epitope.99)), 
                        rep("Non-epitopes Bush 1999", length(dist.non_epitope.99)), 
                        rep("Epitopes Experimental", length(dist.epitope.rename)), 
                        rep("Non-epitopes Experimental", length(dist.non_epitope.rename))))

library(ggplot2)
library(grid)
library(cowplot)

p1 <- ggplot(df[df$id == "Epitopes Bush 1999" | df$id == "Non-epitopes Bush 1999", ], aes(x=distances, fill=factor(id))) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(limits = c(0, 135)) +
  scale_y_continuous(limits = c(0, 0.033)) +
  scale_fill_discrete(name="") +
  xlab('Distances (Angstroms)') +
  ylab('Probability Density') +
  theme(legend.position = c(0.75,0.75))
show(p1)

p2 <- ggplot(df[df$id == "Epitopes Experimental" | df$id == "Non-epitopes Experimental", ], aes(x=distances, fill=factor(id))) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(limits = c(0, 135)) +
  scale_y_continuous(limits = c(0, 0.033)) +
  scale_fill_discrete(name="") +
  xlab('Distances (Angstroms)') +
  ylab('Probability Density') +
  theme(legend.position = c(0.75,0.75))
show(p2)

p3 <- ggplot(df[df$id == "Epitopes Bush 1999" | df$id == "Non-epitopes Bush 1999", ], aes(x=rsa, fill=factor(id))) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 3.75)) +
  scale_fill_discrete(name="") +
  xlab('Relative Solvent Accessibility') +
  ylab('Probability Density') +
  theme(legend.position = c(0.75,0.75))
show(p3)

p4 <- ggplot(df[df$id == "Epitopes Experimental" | df$id == "Non-epitopes Experimental", ], aes(x=rsa, fill=factor(id))) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 3.75)) +
  scale_fill_discrete(name="") +
  xlab('Relative Solvent Accessibility') +
  ylab('Probability Density') +
  theme(legend.position = c(0.75,0.75))
show(p4)

p <- plot_grid(p1, p2, p3, p4, cols=2)
p <- p + draw_plot_label(c("A", "B", "C", "D"), c(0, 1/2, 0, 1/2), c(1, 1, 1/2, 1/2))
ggsave("~/Google Drive/Data/influenza_HA_evolution/distance_plots/distance_distributions.pdf", p, width=10, height=10)