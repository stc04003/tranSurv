## ------------------------------------------------------------------------------------------
## Library and data
## ------------------------------------------------------------------------------------------
library(tranSurv)

data(channing, package = "boot")
chan <- subset(channing, entry < exit)

## ------------------------------------------------------------------------------------------

trReg(Surv(entry, exit, cens) ~ sex, data = chan)
trReg(Surv(entry, exit, cens) ~ sex, data = chan, method = "adjust", control = list(G = 10))


trReg(Surv(entry, exit, cens) ~ sex, data = chan, B = 50)

fit <- trReg(Surv(entry, exit, cens) ~ 1, data = chan)
plot(fit)


(fit <- with(chan, trSurvfit(entry, exit, cens)))
plot(fit)

fit0 <- with(chan, trSurvfit(entry, exit, cens))
str(fit0)



fit <- trReg(Surv(entry, exit, cens) ~ sex, data = chan)
str(fit)

gof(fit)
