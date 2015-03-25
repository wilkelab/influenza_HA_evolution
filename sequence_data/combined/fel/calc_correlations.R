rm(list = ls())

distances <- read.table('distances.dat', sep=',')

rates <- read.table('sites.dat', sep=',', head=T)
df <- read.table('~/Google Drive/Data/influenza_HA_evolution/data_table/numbering_table_unix.csv', head=T, sep=',', stringsAsFactors=F)
map <- df$pdb.4fnk
map[25] <- NA
map[26] <- NA
#dN.dS <- rates$dN
dN.dS <- rates$dN.dS[!is.na(map)]

distances <- distances[-1,]
distances <- distances[-1,]
distances <- distances[,-1]
distances <- distances[,-1]

correlations <- as.vector(sapply(distances, function(x) cor(1/x[x!=0], dN.dS[x!=0])))
p.values <- as.vector(sapply(distances, function(x) cor.test(1/x[x!=0], dN.dS[x!=0])$p.value))
write.table(data.frame(correlations), file= 'all.correlations', row.names=F, col.names=F)
write.table(data.frame(distances[, which(max(dN.dS) == dN.dS)]), file='best_distances.dat', row.names=F, col.names=F)
print(min(correlations))
print(max(correlations))
print(df$Protein[which(max(dN.dS) == dN.dS)])

print(which(min(correlations) == correlations))
print(order(correlations))
