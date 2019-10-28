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


Q <- 3
ti <- seq(min(fit$.data$stop), max(fit$.data$stop), length.out = Q)[-c(1, Q)]
getTL(fit$.data$start, fit$.data$stop, ti, B = 10)


getTL <- function(tt, yy, ti, B) {
    dat <- NULL
    dat$tt <- tt
    dat$yy <- yy
    dat$xt <- dat$yy - dat$tt
    dat <- data.frame(dat)
    nCpt <- length(ti)
    covr <- matrix(NA, nrow = nrow(dat), ncol = nCpt + 1)
    covr[,1] <- with(dat, pmin(yy, ti[1]))
    if (length(ti) > 1) {
        covr[, ncol(covr)] <- pmax(dat$yy - ti[2], 0)
        for (i in 2:nCpt)
            covr[,i] <- pmin(pmax(dat$yy - ti[i - 1], 0), ti[i] - ti[i - 1])
    } else covr[,2] <- pmax(dat$yy - ti[1], 0)
    colnames(covr) <- c(sapply(1:(nCpt + 1), function(x) paste("t", x, sep = "")))                  
    fit <- truncSP::lt(xt ~ covr, data = dat, covar = TRUE, B = B)
    perm <- combn(length(fit$coefficients[-1]), 2)
    con <- matrix(0, ncol = ncol(perm), nrow = length(fit$coefficients[-1]))
    for (i in 1:ncol(perm)) {
        con[perm[1,i], i] <- 1
        con[perm[2,i], i] <- -1
    }     
    dif.vet <- t(con) %*% fit$coefficients[-1]
    dif.cov <- (t(con) %*% fit$covariance[-1, -1]) %*% con
    if (qr(dif.cov)$rank < nrow(dif.cov)) {
        dif.cov <- dif.cov[1:qr(dif.cov)$rank, 1:qr(dif.cov)$rank]
        dif.vet <- dif.vet[1:qr(dif.cov)$rank]
    }
    coefficients <- fit$coefficients
    pval <- 1 - pchisq((t(dif.vet) %*% solve(dif.cov)) %*% dif.vet, length(dif.vet))
    list(coefficients = coefficients, pval = pval)
}
