library(tranSurv)
library(copula)

rho <- iTau(normalCopula(dim = 2), .5) ## convert kendall's tau to pearson's rho

set.seed(1)
u <- rCopula(2000, normalCopula(rho, dim = 2)) ## generate correlated data
u <- qweibull(u, 2, 1) ## assumes t1 and t2 follows some weibull distribution
colnames(u) <- c("t1", "t2")

## check kendall's tau
cor(u, method = "kendall")
kendall(u)

## apply trancation
u <- subset(as.data.frame(u), t1 < t2) 
dim(u) ## ~50% truncation rate because t1 and t2 follows the same distribution
head(u)

attach(u)

condKendall(t1, t2)$PE
wKendall(t1, t2) ## check with condKendall when no weights

## with Perturbation weights
wKendall(t1, t2, rexp(length(t1)))

## SE estimation
set.seed(2)
sd(replicate(500, wKendall(t1, t2, rexp(length(t1))))) ## 0.018
condKendall(t1, t2)$SE ## 0.208

detach(u)

## small size simulation
do <- function() {
    u <- rCopula(2000, normalCopula(rho, dim = 2)) ## generate correlated data
    u <- qweibull(u, 2, 1) ## assumes t1 and t2 follows some weibull distribution
    colnames(u) <- c("t1", "t2")
    u <- subset(as.data.frame(u), t1 < t2)
    out <- with(u, c(sd(replicate(500, wKendall(t1, t2, rexp(length(t1))))) ,
                     condKendall(t1, t2)$SE))
    out
}

set.seed(3)
rowMeans(replicate(100, do())) ## [1] 0.01819335 0.02064378
