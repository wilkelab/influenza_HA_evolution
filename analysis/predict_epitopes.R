rm(list=ls())

linear.eps <- read.table('../epitope_data/linear_eps.counts', head=T)
nonlinear.eps <- read.table('../epitope_data/nonlinear_eps.counts', head=T)
omega <- read.table('../sequence_data/not_structure/combined/fel/rates.out', head=T)
rate4site <- read.table('../sequence_data/not_structure/combined/rate4site/r4s.res')

rsa.monomer <- read.table('../sequence_data/structure/combined/structure/4fnk_monomer.rsa', head=T)
wcn.monomer <- read.table('../sequence_data/structure/combined/structure/4fnk_monomer.contacts', head=T)
rsa.multimer <- read.table('../sequence_data/structure/combined/structure/4fnk_multimer.rsa', head=T)
wcn.multimer <- read.table('../sequence_data/structure/combined/structure/4fnk_multimer.contacts', head=T)
map <- t(read.table('../sequence_data/structure/combined/fel/short.map', sep=','))

#The epitopes already have the first 16 amino acids cut off for numbering
linear.eps.short <- linear.eps$ep_counts[!is.na(as.vector(map[17:566,1]))]
nonlinear.eps.short <- nonlinear.eps$ep_counts[!is.na(as.vector(map[17:566,1]))]
short.omega <- omega[!is.na(map), ]
short.rate4site <- rate4site[!is.na(map), ]

binary.eps <- rep(0, length(nonlinear.eps.short))
binary.eps[nonlinear.eps.short > 0 | linear.eps.short > 0 ] <- 1

all.dat <- data.frame(linear.eps = linear.eps.short, 
                      nonlinear.eps = nonlinear.eps.short, 
                      omega = short.omega$omega, 
                      rate4site = short.rate4site$V3, 
                      rsa.monomer = rsa.monomer$RSA,
                      wcn.monomer = wcn.monomer$wcn,
                      rsa.multimer = rsa.multimer$RSA,
                      wcn.multimer = wcn.multimer$wcn,
                      rsa.diff = rsa.monomer$RSA - rsa.multimer$RSA,
                      binary.eps = binary.eps
                      )

training.dat <- data.frame(eps = factor(binary.eps),
                           rate4site=short.rate4site$V3,
                           wcn.monomer = wcn.monomer$wcn,
                           wcn.multimer = wcn.multimer$wcn
                           )

#fit <- glm(binary.eps ~ omega + wcn + rsa + rate4site, data=all.dat, family = binomial(link=probit))

#library(fmsb)

# odds.dat <- c(sum(all.dat$nonlinear.eps > 0 & all.dat$omega > 1),  
#               sum(!(all.dat$nonlinear.eps > 0) & all.dat$omega > 1), 
#               sum(all.dat$nonlinear.eps > 0 & !(all.dat$omega > 1)),
#               sum(!(all.dat$nonlinear.eps > 0) & !(all.dat$omega > 1)))
# 
# print(oddsratio(odds.dat[1], odds.dat[2], odds.dat[3], odds.dat[4]))

library(caret)

preProcValues <- preProcess(x = training.dat[, 2:length(training.dat)], method = c("center", "scale", "knnImpute"))
predictions <- predict(preProcValues, training.dat[, 2:length(training.dat)])

train_control <- trainControl(method="repeatedcv", 
                              repeats=10)

training.dat <- cbind(eps=training.dat$eps, predictions)

model <- train(eps~., data=training.dat, 
               trControl=train_control, 
               method="knn"
               )

predictions <- predict(model, training.dat[, 2:length(training.dat)])
confusion.dat <- confusionMatrix(predictions, training.dat$eps)
print(confusion.dat)
