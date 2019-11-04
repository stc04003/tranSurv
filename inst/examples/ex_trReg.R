data(channing, package = "boot")
chan <- subset(channing, entry < exit)
trReg(Surv(entry, exit, cens) ~ sex, data = chan)
trReg(Surv(entry, exit, cens) ~ sex, data = chan, method = "adjust", control = list(G = 10))


library(tranSurv)

fit <- trReg(Surv(entry, exit, cens) ~ sex, data = chan)
gg <- gof(fit, B = 10)


debug(trReg)
debug(trFit.adjust)
trFit.adjust(DF, engine, stdErr)


lapply(split(x$.data, cut(x$.data$stop, ti)), function(d) {
    tmp <- trReg(Surv(start, stop, status) ~ ., data = d, method = x$method)
    tmp$.data$trans <- with(tmp$.data, x$tFun(stop, start, tmp$a))
    tmp$.data$a <- tmp$a
    return(tmp$.data)
})

dat1 <- split(x$.data, cut(x$.data$stop, ti))[[1]]
trReg(Surv(start, stop, status) ~ dat1[,x$vNames], data = dat1,
      control = list(sc = list(time = sc$time, surv = sc$surv)))
