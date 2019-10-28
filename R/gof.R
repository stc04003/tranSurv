#' Goodness of fit
#'
#' @param x an object of class \code{trSurvfit} returned by the \code{trSurvfit()} or the \code{trReg()} function.
#' @param B a numerical value specifies the bootstrap size.
#' @param Q a numerical value specifies the number of cupoints.
#' 
#' @export
#' @example inst/examples/ex_gof.R
#' 
gof.trReg <- function(x, B = 200,Q) {
    B <- max(x$B, B)
    ti <- seq(min(x$.data$stop), max(x$.data$stop), length.out = Q)[-c(1, Q)]
    ltfit <- getTL(x$.data$start, x$.data$stop, ti, B)

    ## print p-values
    ## combine all subsets
}

#' @importFrom truncSP lt
#'
#' @param tt is the truncation time
#' @param yy is the observed survival times (events only)
#' @param ti is the endpoints of grids
#' @param B is the bootstrap size
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
    fit <- lt(xt ~ covr, data = dat, covar = TRUE, B = B)
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
