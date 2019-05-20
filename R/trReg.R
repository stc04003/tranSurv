## Functions for different methods

#' 1. Find quasi-independent latent truncation time by maximizing conditional Kendall's tau based on the uncensored observations only
#' 2. Use the latent truncation time in the Cox model as a truncation time
#'
#' @noRd
#' @keywords internal
trFit.kendall <- function(DF, engine) {
    trun <- DF$start
    obs <- DF$stop
    delta <- DF$status
    ## Finding the latent truncation time
    lower <- ifelse(engine@lower == -Inf,  -.Machine$integer.max, engine@lower)
    upper <- ifelse(engine@upper == Inf, .Machine$integer.max, engine@upper)
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    trun1 <- trun[order(obs)][delta[order(obs)] == 1]
    obs1 <- sort(obs[delta == 1])
    delta1 <- delta[delta == 1]
    ## optimize getA in many grids
    grids <- seq(lower + 1e-5, upper, length.out = engine@G)
    tmp <- sapply(1:(engine@G - 1), function(y)
        optimize(f = function(x) abs(getA(x, trun1, obs1, delta1, sc = sc, FUN = engine@tFun)$PE),
                 tol = engine@tol, interval = c(grids[y], grids[y + 1])))
    a <- as.numeric(tmp[1, which.min(tmp[2,])])
    ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
    wgtX <- approx(sc$time, sc$surv, obs1, "constant", yleft = 1, yright = min(sc$surv))$y
    out <- coef(summary(coxph(Surv(ta, obs1, delta1) ~ DF[delta == 1,-(1:3)], weights = 1 / wgtX)))    
    ## trun0 <- ifelse(delta == 1, ta[match(obs, obs1)], trun)
    ## out <- coef(summary(coxph(Surv(trun0, obs, delta) ~ DF[,-(1:3)], weights = 1 / wgtX)))
    rownames(out) <- names(DF)[-(1:3)]
    return(out)
}

trFit.adjust <- function(DF, engine) {
    trun <- DF$start
    obs <- DF$stop
    delta <- DF$status
    lower <- ifelse(engine@lower == -Inf,  -.Machine$integer.max, engine@lower)
    upper <- ifelse(engine@upper == Inf, .Machine$integer.max, engine@upper)
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    trun1 <- trun[order(obs)][delta[order(obs)] == 1]
    obs1 <- sort(obs[delta == 1])
    delta1 <- delta[delta == 1]
    wgtX <- approx(sc$time, sc$surv, obs1, "constant", yleft = 1, yright = min(sc$surv))$y
    coxAj <- function(a) {
        ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
        if (engine@Q > 0)
            cov <- model.matrix( ~ cut(ta, breaks = quantile(ta, 0:(1 + engine@Q) / (1 + engine@Q)),
                                       include.lowest = TRUE) - 1)
        else cov <- ta
        min(sum(coef(coxph(Surv(ta, obs1, delta1) ~ DF[delta == 1,-(1:3)] + cov, weights = 1 / wgtX))[-(1:(NCOL(DF) - 3))]^2, na.rm = TRUE), 1e4)
    }
    grids <- seq(lower + 1e-5, upper, length.out = engine@G)
    tmp <- sapply(1:(engine@G - 1), function(y)
        optimize(f = function(x) suppressWarnings(coxAj(x)), interval = c(grids[y], grids[y + 1])))
    a <- as.numeric(tmp[1, which.min(tmp[2,])])
    ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
    if (engine@Q > 0)
        cov <- model.matrix( ~ cut(ta, breaks = quantile(ta, 0:(1 + engine@Q) / (1 + engine@Q)),
                                   include.lowest = TRUE) - 1)
    else cov <- ta
    out <- coef(summary(coxph(Surv(ta, obs1, delta1) ~ DF[delta == 1,-(1:3)] + cov, weights = 1 / wgtX)))[1:(NCOL(DF) - 3),,drop = FALSE]
    rownames(out) <- names(DF)[-(1:3)]
    return(out)    
}

## Class definition
setClass("Engine",
         representation(tol = "numeric", lower = "numeric", upper = "numeric", G = "numeric", Q = "numeric", tFun = "function"),
         prototype(tol = 1e-2, lower = -1, upper = 50, G = 30, Q = 0),
         contains= "VIRTUAL")
setlass("kendal", contains = Engine)
setlass("adjust", contains = Engine)
                        
## Method Dispatch
setGeneric("trFit", function(DF, engine) {standardGeneric("trFit")})

setMethod("trFit", signature(engine = "kendall"), trFit.kendall)
setMethod("trFit", signature(engine = "adjust"), trFit.adjust)

#' Fitting regression model via structural transformation model
#'
#' \code{trReg} fits transformation model under dependent truncation and independent censoring via a structural transformation model.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#' The \code{response} is assumed to be a \code{survival::Surv} object with both left truncation and right censoring.
#' See \code{?survival::Surv} for more details.
#' @param data  an optional data.frame in which to interpret the variables occurring
#'     in the \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param tFun a character string specifying the transformation function or a user specified function indicating the relationship between X, T, and a.
#' When \code{tFun} is a character, the following are permitted:
#' \describe{
#'   \item{linear}{linear transformation structure}
#'   \item{log}{log-linear transformation structure}
#'   \item{exp}{exponential transformation structure}
#' }
#' @param method a character string specifying the underlying model. See \bold{Details}.
#'
#' @importFrom survival is.Surv coxph
#' @export
#' 
trReg <- function(formula, data, subset, tFun = "linear",
                  method = c("kendall", "adjust")) {
    if (class(tFun) == "character") {
        if (tFun == "linear") FUN <- function(X, T, a) (T + a * X) / (1 + a)
        if (tFun == "log") FUN <- function(X, T, a) exp((log(replace(T, 0, 1)) + a * log(X)) / (1 + a))
        if (tFun == "log2") FUN <- function(X, T, a) exp((1 + a) * log(replace(T, 0, 1)) - a * log(X))
        if (tFun == "exp") FUN <- function(X, T, a) log((exp(T) + a * exp(X)) / (1 + a))
    } else {
        FUN <- match.fun(tFun)
    }
    method <- match.arg(method)
    Call <- match.call()
    if (missing(data)) {
        res <- eval(formula[[2]], parent.frame())
        cov <- model.matrix(formula, parent.frame())
    } else {
        res <- eval(formula[[2]], data)
        cov <- model.matrix(formula, data)
    }   
    if (!is.Surv(res)) stop("Response must be a Surv resect")
    if (!match("start", attr(res, "dimnames")[[2]])) stop("Missing left-truncation time")
    formula[[2]] <- NULL
    if (formula == ~1) DF <- as.data.frame(unclass(res))
    else DF <- as.data.frame(cbind(res, cov)) ## First 3 columns reserved to `start`, `stop`, `status`
    DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    out <- trFit(DF, engine)
    out
}
