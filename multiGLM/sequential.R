## run sequential testing
## load data
source("sequential.dat.R")

## define SNP and trait matrix
snp <- sequential.dat[,1]
y <-  sequential.dat[,2:17]
n.traits <- ncol(y)
glm.family <- rep("gaussian", n.traits)

# make sure pleio package is intalled, available in CRAN
library(pleio)

## first, fit to get summary statistics
fit <- pleio.glm.fit(y, snp, glm.family)

## second, perform sequential testing
stat.seq <- pleio.glm.sequential(fit, pval.threshold=.05)
stat.seq

## third, look at each step

stat0 <- pleio.glm.test(fit, count.nonzero.coef = 0)
stat0$pval

stat1 <- pleio.glm.test(fit, count.nonzero.coef = 1)
stat1$pval

stat2 <- pleio.glm.test(fit, count.nonzero.coef = 2)
stat2$pval

stat3 <- pleio.glm.test(fit, count.nonzero.coef = 3)
stat3$pval

stat4 <- pleio.glm.test(fit, count.nonzero.coef = 4)
stat4$pval
