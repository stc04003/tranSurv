## Functions for different methods

#' 1. Find quasi-independent latent truncation time by maximizing conditional Kendall's tau based on the uncensored observations only
#' 2. Use the latent truncation time in the Cox model as a truncation time
#'
#' @noRd
#' @keywords internal
trFit.kendall <- function(DF, engine, stdErr) {
    out <- NULL
    trun <- DF$start
    obs <- DF$stop
    delta <- DF$status
    ## Finding the latent truncation time
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    trun1 <- trun[delta == 1] ## trun[order(obs)][delta[order(obs)] == 1]
    obs1 <- obs[delta == 1]
    delta1 <- delta[delta == 1]
    ## optimize getA in many grids
    grids <- seq(engine@lower + 1e-5, engine@upper, length.out = engine@G)
    tmp <- sapply(1:(engine@G - 1), function(y)
        optimize(f = function(x) abs(getA(x, trun1, obs1, delta1, sc = sc, FUN = engine@tFun)$PE),
                 tol = engine@tol, interval = c(grids[y], grids[y + 1])))
    a <- as.numeric(tmp[1, which.min(tmp[2,])])
    ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
    ## guard against infity weights
    if (min(sc$surv) == 0) 
        sc$surv <- ifelse(sc$surv == min(sc$surv), sort(unique(sc$surv))[2], sc$surv)
    wgtX <- approx(sc$time, sc$surv, obs1, "constant", yleft = 1, yright = min(sc$surv))$y
    out$PE <- coef(summary(coxph(Surv(ta, obs1, delta1) ~ as.matrix(DF[delta == 1,-(1:3)]),
                                 weights = 1 / wgtX)))
    ## trun0 <- ifelse(delta == 1, ta[match(obs, obs1)], trun)
    ## out <- coef(summary(coxph(Surv(trun0, obs, delta) ~ DF[,-(1:3)], weights = 1 / wgtX)))
    rownames(out$PE) <- names(DF)[-(1:3)]
    out$varNames <- names(DF)[-(1:3)]
    out$SE <- NA
    out$a <- a
    return(out)
}

#' @noRd
#' @keywords internal
trFit.adjust <- function(DF, engine, stdErr) {
    out <- NULL
    trun <- DF$start
    obs <- DF$stop
    delta <- DF$status
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    trun1 <- trun[delta == 1]
    obs1 <- obs[delta == 1]
    delta1 <- delta[delta == 1]
    if (min(sc$surv) == 0) 
        sc$surv <- ifelse(sc$surv == min(sc$surv), sort(unique(sc$surv))[2], sc$surv)
    wgtX <- approx(sc$time, sc$surv, obs1, "constant", yleft = 1, yright = min(sc$surv))$y
    coxAj <- function(a) {
        ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
        if (engine@Q > 0)
            cov <- model.matrix( ~ cut(ta, breaks = quantile(ta, 0:(1 + engine@Q) / (1 + engine@Q)),
                                       include.lowest = TRUE) - 1)
        else cov <- ta
        min(sum(coef(coxph(Surv(ta, obs1, delta1) ~ as.matrix(DF[delta == 1,-(1:3)]) + cov,
                           weights = 1 / wgtX))[-(1:(NCOL(DF) - 3))]^2, na.rm = TRUE), 1e4)
    }
    grids <- seq(engine@lower + 1e-5, engine@upper, length.out = engine@G)
    tmp <- sapply(1:(engine@G - 1), function(y)
        optimize(f = function(x) suppressWarnings(coxAj(x)), interval = c(grids[y], grids[y + 1])))
    a <- as.numeric(tmp[1, which.min(tmp[2,])])
    ta <- mapply(engine@tFun, X = obs1, T = trun1, a = a)
    if (engine@Q > 0)
        cov <- model.matrix( ~ cut(ta, breaks = quantile(ta, 0:(1 + engine@Q) / (1 + engine@Q)),
                                   include.lowest = TRUE) - 1)
    else cov <- ta
    out$PE <- coef(summary(coxph(Surv(ta, obs1, delta1) ~ as.matrix(DF[delta == 1,-(1:3)]) + cov,
                              weights = 1 / wgtX)))[1:(NCOL(DF) - 3),,drop = FALSE]
    rownames(out$PE) <- names(DF)[-(1:3)]
    out$SE <- NA
    out$varNames <- names(DF)[-(1:3)]
    out$a <- a
    return(out)    
}

#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#' @noRd
#' @keywords internal
trFit.boot <- function(DF, engine, stdErr) {
    trun <- DF$start
    obs <- DF$stop
    delta <- DF$status
    out <- trFit(DF, engine, NULL)
    win <- (engine@upper - engine@lower) / engine@G
    engine@lower <- max(engine@lower, out$a - win)
    engine@upper <- min(engine@upper, out$a + win)
    engine@G <- max(5, round(engine@G / 10))
    if (stdErr@parallel) {
        cl <- makeCluster(stdErr@parCl)
        clusterExport(cl = cl,
                      varlist = c("DF", "engine"), envir = environment())
        out$SE <- parSapply(cl, 1:stdErr@B, function(x)
            trFit(DF[sample(1:NROW(DF), NROW(DF), TRUE),], engine, NULL)$PE[,1])
        stopCluster(cl)
    } else out$SE <- replicate(stdErr@B,
                               trFit(DF[sample(1:NROW(DF), NROW(DF), TRUE),], engine, NULL)$PE[,1])
    if (nrow(out$PE) > 1) out$SE <- apply(out$SE, 1, sd)
    else out$SE <- sd(out$SE)
    out
}

#' Class definition
#' @noRd
#' @keywords internal
setClass("Engine",
         representation(tol = "numeric", lower = "numeric", upper = "numeric",
                        G = "numeric", Q = "numeric", tFun = "function"),
         prototype(tol = 1e-2, lower = -1, upper = 20, G = 50, Q = 0),
         contains= "VIRTUAL")
setClass("kendall", contains = "Engine")
setClass("adjust", contains = "Engine")

setClass("stdErr",
         representation(B = "numeric", parallel = "logical", parCl = "numeric"),
         prototype(B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2),
         contains = "VIRTUAL")
setClass("bootstrap", contains = "stdErr")

#' Method Dispatch
#' @noRd
#' @keywords internal
setGeneric("trFit", function(DF, engine, stdErr) {standardGeneric("trFit")})

setMethod("trFit", signature(engine = "kendall", stdErr = "NULL"), trFit.kendall)
setMethod("trFit", signature(engine = "adjust", stdErr = "NULL"), trFit.adjust)
setMethod("trFit", signature(engine = "kendall", stdErr = "bootstrap"), trFit.boot)
setMethod("trFit", signature(engine = "adjust", stdErr = "bootstrap"), trFit.boot)

#' Fitting regression model via structural transformation model
#'
#' \code{trReg} fits transformation model under dependent truncation and independent censoring via a structural transformation model.
#'
#' The main assumption on the structural transformation model is that it assumes there is a latent, quasi-independent truncation time
#' that is associated with the observed dependent truncation time, the event time, and an unknown dependence parameter
#' through a specified funciton.
#' The structure of the transformation model is of the form:
#' \deqn{h(U) = (1 + a)^{-1} \times (h(T) + ah(X)),} where \eqn{T} is the truncation time, \eqn{X} is the observed failure time,
#' \eqn{U} is the transformed truncation time that is quasi-independent from \eqn{X} and \eqn{h(\cdot)} is a monotonic transformation function.
#' The condition, \eqn{T < X}, is assumed to be satisfied.
#' The quasi-independent truncation time, \eqn{U}, is obtained by inverting the test for quasi-independence by one of the following methods:
#' \describe{
#'   \item{\code{method = "kendall"}}{ by minimizing the absolute value of the restricted inverse probability weighted Kendall's tau or maximize the corresponding \eqn{p}-value.
#' This is the same procedure used in the \code{trSUrvfit()} function.}
#'   \item{\code{method = "adjust"}}{ includes a function of latent truncation time, \eqn{U}, as a covariate.
#'   A piece-wise function is constructed based on (\eqn{Q + 1}) indicator functions on whether \eqn{U} falls in the \eqn{Q}th and the (\eqn{Q+1})th percentile,
#'   where \eqn{Q} is the number of cutpoints used. See \code{control} for details. 
#'   The transformation parameter, \eqn{a}, is then chosen to minimize the significance of the coefficient parameter.}} 
#' 
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#' The \code{response} is assumed to be a \code{survival::Surv} object with both left truncation and right censoring.
#' When there is no covariates, e.g., when the right hand side of the formula is \code{~ 1}, the \code{trReg()} function returns a \code{trSurvfit} object.
#' See \code{?survival::Surv} for more details.
#' @param data  an optional data frame in which to interpret the variables occurring
#'     in the \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param tFun a character string specifying the transformation function or a user specified function indicating the relationship between \eqn{X}, \eqn{T}, and \eqn{a}.
#' When \code{tFun} is a character, the following are permitted:
#' \describe{
#'   \item{linear}{linear transformation structure,}
#'   \item{log}{log-linear transformation structure,}
#'   \item{exp}{exponential transformation structure.}
#' }
#' @param method a character string specifying how the transformation parameter is estimated. The available options are \code{"kendall"} and \code{"adjust"}. See \bold{Details}.
#' @param B a numerical value specifies the bootstrap size for estimating the standard error.
#' When \code{B = 0} (default), the bootstrap standard errors will not be computed.
#' @param control a list of control parameters. The following arguments are allowed:
#' \describe{
#'   \item{\code{lower}}{The lower bound to search for the transformation parameter; default at -1.}
#'   \item{\code{upper}}{The upper bound to search for the transformation parameter; default at 20.}
#'   \item{\code{tol}}{The tolerance used in the search for the transformation parameter; default at 0.01.}
#'   \item{\code{G}}{The number of grids used in the search for the transformation parameter; default is 50.
#' A smaller \code{G} could results in faster search, but might be inaccurate.}
#'   \item{\code{Q}}{The number of cutpoints for the truncation time used when \code{method = "adjust"}. The default is 0.}
#'   \item{\code{parallel}}{an logical value indicating whether parallel computation will be applied when \code{B > 0}.}
#'   \item{\code{parCl}}{an integer value specifying the number of CPU cores to be used when \code{parallel = TRUE}.
#' The default value is half the CPU cores on the current host.}
#' }
#'
#' 
#' @importFrom survival is.Surv coxph
#' @importFrom methods getClass
#' @seealso \code{\link{trSurvfit}}
#' 
#' @export
#' @example inst/examples/ex_trReg.R
trReg <- function(formula, data, subset, tFun = "linear",
                  method = c("kendall", "adjust"),
                  B = 0, control = list()) {
    method <- match.arg(method)
    Call <- match.call()
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class = method), engine.control))
    stdErr.control <- control[names(control) %in% names(attr(getClass("bootstrap"), "slots"))]
    stdErr <- do.call("new", c(list(Class = "bootstrap"), stdErr.control))
    stdErr@B <- B
    if (B == 0) class(stdErr)[[1]] <- "NULL"
    if (class(tFun) == "character") {
        if (tFun == "linear") engine@tFun <- function(X, T, a) (T + a * X) / (1 + a)
        if (tFun == "log") engine@tFun <- function(X, T, a) exp((log(replace(T, 0, 1)) + a * log(X)) / (1 + a))
        if (tFun == "log2") engine@tFun <- function(X, T, a) exp((1 + a) * log(replace(T, 0, 1)) - a * log(X))
        if (tFun == "exp") engine@tFun <- function(X, T, a) log((exp(T) + a * exp(X)) / (1 + a))
    } else {
        engine@tFun <- match.fun(tFun)
    }
    if (missing(data)) {
        resp <- eval(formula[[2]], parent.frame())
        covM <- model.matrix(formula, parent.frame())
    } else {
        resp <- eval(formula[[2]], data)
        covM <- model.matrix(formula, data)
    }
    vNames <- setdiff(colnames(covM), "(Intercept)") # attr(terms(formula), "term.labels")
    if (!is.Surv(resp)) stop("Response must be a Surv resect")
    if (!match("start", attr(resp, "dimnames")[[2]])) stop("Missing left-truncation time")
    engine@lower <- ifelse(engine@lower == -Inf,  -.Machine$integer.max, engine@lower)
    engine@upper <- ifelse(engine@upper == Inf, .Machine$integer.max, engine@upper)
    formula[[2]] <- NULL
    if (formula == ~1) {
        DF <- as.data.frame(unclass(resp))
        out <- trSurvfit(DF$start, DF$stop, DF$status, trans = tFun, plots = FALSE,
                         control = trSurv.control(lower = engine@lower, upper = engine@upper))
        class(out) <- "trSurvfit"
    } else {
        DF <- as.data.frame(cbind(resp, covM)) ## First 3 columns reserved to `start`, `stop`, `status`
        DF <- DF[,which(colnames(DF) != "(Intercept)")]
        out <- trFit(DF, engine, stdErr)
        class(out) <- "trReg"
    }
    out$Call <- Call
    out$B <- B
    out$Q <- engine@Q
    out$tFun <- engine@tFun
    out$vNames <- vNames
    out$method <- method
    out$.data <- DF 
    out
}
