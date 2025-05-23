% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trReg.R
\name{trReg}
\alias{trReg}
\title{Fitting regression model via structural transformation model}
\usage{
trReg(
  formula,
  data,
  subset,
  tFun = "linear",
  method = c("kendall", "adjust"),
  B = 0,
  control = list()
)
}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}.
The \code{response} is assumed to be a \code{survival::Surv} object with both left truncation and right censoring.
When there is no covariates, e.g., when the right hand side of the formula is \code{~ 1}, the \code{trReg()} function returns a \code{trSurvfit} object.
See \code{?survival::Surv} for more details.}

\item{data}{an optional data frame in which to interpret the variables occurring
in the \code{formula}.}

\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}

\item{tFun}{a character string specifying the transformation function or a user specified function indicating the relationship between \eqn{X}, \eqn{T}, and \eqn{a}.
When \code{tFun} is a character, the following are permitted:
\describe{
  \item{linear}{linear transformation structure,}
  \item{log}{log-linear transformation structure,}
  \item{exp}{exponential transformation structure.}
}}

\item{method}{a character string specifying how the transformation parameter is estimated. The available options are \code{"kendall"} and \code{"adjust"}. See \bold{Details}.}

\item{B}{a numerical value specifies the bootstrap size for estimating the standard error.
When \code{B = 0} (default), the bootstrap standard errors will not be computed.}

\item{control}{a list of control parameters. The following arguments are allowed:
\describe{
  \item{\code{lower}}{The lower bound to search for the transformation parameter; default at -1.}
  \item{\code{upper}}{The upper bound to search for the transformation parameter; default at 20.}
  \item{\code{tol}}{The tolerance used in the search for the transformation parameter; default at 0.01.}
  \item{\code{G}}{The number of grids used in the search for the transformation parameter; default is 50.
A smaller \code{G} could results in faster search, but might be inaccurate.}
  \item{\code{Q}}{The number of cutpoints for the truncation time used when \code{method = "adjust"}. The default is 0.}
  \item{\code{P}}{The number of breakpoints to divide the event times into equally spaced segments.
When \code{P > 1}, the latent truncation time, \eqn{T'(a)} will be computed in each subset.
The transformation model is then applied to the aggregated data.} 
  \item{\code{a}}{The transformation parameter. When this is specified, the transformation model is applied based on the specified \code{a}.
                  When this is not specified, an optimized \code{a} will be determined by optimization one of the quasi-independence measure. See \bold{Details}.}
  \item{\code{parallel}}{an logical value indicating whether parallel computation will be applied when \code{B > 0}.}
  \item{\code{parCl}}{an integer value specifying the number of CPU cores to be used when \code{parallel = TRUE}.
The default value is half the CPU cores on the current host.}
}}
}
\value{
A \code{trReg} object containing the following components:
\describe{
  \item{\code{PE}}{A named numeric matrix of point estimates and related statistics (e.g., coefficient, exponentiated coefficient, standard error, z-score, and p-value).}
  \item{\code{varNames}}{Character string giving the name(s) of the covariates.}
  \item{\code{SE}}{A numeric vector contains the bootstrap standard error.}
  \item{\code{a}}{Estimated transformation parameter.}
  \item{\code{Call}}{The matched call to the fitting function.}
  \item{\code{B}, \code{Q}, \code{P}}{Model parameters; \code{B} is the bootstrap sapmle, \code{Q} is the number of cutpoints, and \code{P} is the number of break points. See \bold{Details}.}
  \item{\code{tFun}}{A function defining the transformation model.}
  \item{\code{vNames}}{Character vector of covariate names.}
  \item{\code{method}}{Character string specifying the estimation method (e.g., \code{"kendall"} or \code{"adjust"}).}
  \item{\code{.data}}{A data frame used in fitting.}
}
}
\description{
\code{trReg} fits transformation model under dependent truncation and independent censoring via a structural transformation model.
}
\details{
The main assumption on the structural transformation model is that it assumes there is a latent, quasi-independent truncation time
that is associated with the observed dependent truncation time, the event time, and an unknown dependence parameter
through a specified function.
The structure of the transformation model is of the form:
\deqn{h(U) = (1 + a)^{-1} \times (h(T) + ah(X)),} where \eqn{T} is the truncation time, \eqn{X} is the observed failure time,
\eqn{U} is the transformed truncation time that is quasi-independent from \eqn{X} and \eqn{h(\cdot)} is a monotonic transformation function.
The condition, \eqn{T < X}, is assumed to be satisfied.
The quasi-independent truncation time, \eqn{U}, is obtained by inverting the test for quasi-independence by one of the following methods:
\describe{
  \item{\code{method = "kendall"}}{ by minimizing the absolute value of the restricted inverse probability weighted Kendall's tau or maximize the corresponding \eqn{p}-value.
This is the same procedure used in the \code{trSUrvfit()} function.}
  \item{\code{method = "adjust"}}{ includes a function of latent truncation time, \eqn{U}, as a covariate.
  A piece-wise function is constructed based on (\eqn{Q + 1}) indicator functions on whether \eqn{U} falls in the \eqn{Q}th and the (\eqn{Q+1})th percentile,
  where \eqn{Q} is the number of cutpoints used. See \code{control} for details. 
  The transformation parameter, \eqn{a}, is then chosen to minimize the significance of the coefficient parameter.}}
}
\examples{
data(channing, package = "boot")
chan <- subset(channing, entry < exit)
trReg(Surv(entry, exit, cens) ~ sex, data = chan)

}
\seealso{
\code{\link{trSurvfit}}
}
