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
#' @param tFun a character string specifying the transformation function or a user specified function.
#' When \code{tFun} is a character, the following are permitted:
#' \describe{
#'   \item{linear}{linear transformation structure; equivalent to \code{function(x) x}}
#'   \item{log}{log-linear transformation structure; equivalent to \code{function(x) log(x)}}
#'   \item{exp}{exponential transformation structure; equivalent to \code{function(x) exp(x)}}
#' }
#' @param method a character string specifying the underlying model. See \bold{Details}.
#'
#' @importFrom survival is.Surv
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
    DF
}
