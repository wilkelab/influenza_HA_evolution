rm(list = ls())
setwd('~/Google Drive/Data/influenza_HA_evolution/manuscript/')
df <- read.table('numbering_table.csv', head=T, sep=',')
df <- df[!is.na(df$pdb.4fnk) & !is.na(df$distance.to.224), ]

library(boot)

fit <- glm(FEL.dN.dS ~ I(1/distance.to.224), data = df)
df <- cbind(df, predicted.w=predict(fit, df))
r.value <- cor(df$FEL.dN.dS, df$predicted.w, use="complete.obs")
rse <- sqrt(deviance(fit)/df.residual(fit))
cv.fit <- cv.glm(data=df, glmfit=fit, K=10)
cv.rse <- cv.fit$delta

print(r.value^2)
print(rse)
print(sqrt(cv.rse[1]))
