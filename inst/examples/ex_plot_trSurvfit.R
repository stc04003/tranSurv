data(channing, package = "boot")
chan <- subset(channing, entry < exit)

plot(trReg(Surv(entry, exit, cens) ~ 1, data = chan))
plot(with(chan, trSurvfit(entry, exit, cens)))
plot(trReg(Surv(entry, exit, cens) ~ sex, data = chan))
plot(trReg(Surv(entry, exit, cens) ~ sex, data = chan, method = "adjust", control = list(G = 10)))
