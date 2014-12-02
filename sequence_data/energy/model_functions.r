remove(list=ls())

setwd('~/Google Drive/Data/influenza_HA_evolution/sequence_data/energy/')

require(MASS)
# Stability Model
# Derived following Bloom (2009) assumptions
# Assuming an unbiased aaJC amino-acid mutation model
# It is a Monte Carlo Process 
# Self-Consistent with the equilibrium frequencies
predict.rate.stability_model = function( ddG.row, scale = 1 )
{
    pf = function(x) min(1,x)

    x = unlist(ddG.row[2:21])
    p = exp(-scale*(x-min(x)))
    p.norm = p/sum(p) # equilibrium distribution for the previous qmat
    rate = 2*sum((rank(-p.norm)-1)*p.norm)
    rate
}

d <- read.table('4FNK_A_foldx_ddG.txt', head=T)

ddg <- data.frame(ddg=apply(d, 1, predict.rate.stability_model))

write.table(ddg, file = 'ddG.txt', row.names = F, col.names = F)