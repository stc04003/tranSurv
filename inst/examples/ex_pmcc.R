## Generate simulated data from transformation model
data(channing, package = "boot")
chan <- subset(channing, sex == "Male" & entry < exit & cens == 1)
with(chan, pmcc(entry, exit)) ## cannot handle right censored data
