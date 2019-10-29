library(tranSurv)

## Generate simulated data from a transformation model
datgen <- function(n) {
    a <- -0.3
    X <- rweibull(n, 2, 4) ## failure times
    U <- rweibull(n, 2, 1) ## latent truncation time
    T <- (1 + a) * U - a * X ## apply transformation
    C <- 10 ## censoring
    dat <- data.frame(trun = T, obs = pmin(X, C), delta = 1 * (X <= C))
    return(subset(dat, trun <= obs))
}

set.seed(123)
dat <- datgen(300)
fit <- with(dat, trSurvfit(trun, obs, delta))
fit

## Checking the transformation parameter
fit$byTau$par
fit$byTau$obj
with(dat, condKendall(trun, obs, delta, method = "IPW2", a = fit$byTau$par[1]))$PE

fit$byP$par
fit$byP$obj
with(dat, condKendall(trun, obs, delta, method = "IPW2", a = fit$byP$par[1]))$p.value



for (i in 1:1000) {
    set.seed(i)
    dat <- datgen(300)
    fit <- with(dat, trSurvfit(trun, obs, delta))
}

## ---------------------------------------------------------------------------------
## Debugging
## ---------------------------------------------------------------------------------

set.seed(123)
dat <- datgen(300)
debug(trSurvfit)

with(dat, trSurvfit(trun, obs, delta))




uniroot.all0 <- function (f, interval, lower= min(interval),
                         upper= max(interval), ... ) {
    n = 100
    xseq <- seq(lower, upper, len = n + 1)
    mod  <- f(xseq,...)
    Equi <- xseq[which(mod==0)]
    ss   <- mod[1:n] * mod[2:(n+1)]
    ii   <- which(ss<0)
    ## print(c(mod[ii], mod[ii + 1]))
    if (length(ii) > 0) 
    for (i in ii) Equi <- c(Equi,uniroot(f,lower=xseq[i],upper=xseq[i+1],...)$root)
    else Equi <- NULL
    return(Equi)
}


fun <- function (x) cos(2*x)^3

curve(fun(x), 0, 10,main = "uniroot.all")

uniroot.all0(fun, c(0, 10))
uniroot.all0(fun, c(.8, 1.8))

uniroot.all(fun, c(0, 10))
uniroot.all(fun, c(.8, 1.8))



## a difficult case...
f <- function (x) 1/cos(1+x^2)
AA <- uniroot.all(f, c(-5, 5))
curve(f(x), -5, 5, n = 500, main = "uniroot.all")
points(AA, rep(0, length(AA)), col = "red", pch = 16)












f1 <- function(x) sapply(ii, function(x) uniroot(f = fun, lower = xseq[x], upper = xseq[x + 1], ...)$root)
f2 <- function(x) {
    for (i in ii) Equi <- c(Equi,uniroot(f,lower=xseq[i],upper=xseq[i+1],...)$root)
    Equi
}

identical(f1(), f2())
library(microbenchmark)
microbenchmark(f1(), f2(), times = 1e3)
