data(channing, package = "boot")
chan <- subset(channing, entry < exit)
trReg(Surv(entry, exit, cens) ~ sex, data = chan)
trReg(Surv(entry, exit, cens) ~ sex, data = chan, method = "adjust", control = list(G = 10))

