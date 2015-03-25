rm(list = ls())

distances <- read.table('1HTM_distance_matrix.dat', sep=',')

dat <- read.table('~/Google Drive/Data/influenza_HA_evolution/data_table/numbering_table_unix.csv', head=T, sep=',')
map <- dat$pdb.4fnk
map[25] <- NA
map[26] <- NA
dN.dS <- dat$FEL.dN.dS[!is.na(dat$RSA.Postfusion)]

correlations <- as.vector(sapply(distances, function(x) cor(1/x[x!=0], dN.dS[x!=0])))
p.values <- as.vector(sapply(distances, function(x) cor.test(1/x[x!=0], dN.dS[x!=0])$p.value))
p.values <- p.adjust(p.values, method='fdr')
write.table(data.frame(correlations), file= 'internals.correlations', row.names=F, col.names=F)
print(min(correlations))
print(max(correlations))
print(which(max(dN.dS) == dN.dS))

print(which(min(correlations) == correlations))
print(order(correlations))
print(p.values)
