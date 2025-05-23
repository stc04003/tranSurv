% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kendalls.R
\name{cKendall}
\alias{cKendall}
\title{Conditional Kendall's tau}
\usage{
cKendall(
  trun,
  obs,
  delta = NULL,
  method = "MB",
  weights = NULL,
  a = 0,
  trans = "linear",
  ...
)
}
\arguments{
\item{trun}{left truncation time satisfying \code{trun} <= \code{obs}.}

\item{obs}{observed failure time, must be the same length as \code{trun}, might be right-censored.}

\item{delta}{an optional 0-1 vector of censoring indicator (0 = censored, 1 = event) for \code{obs}.
If this vector is not specified, \code{cKendall} assumes no censoring and all observed failure time
denote events.}

\item{method}{a character string specifying the different version of conditional Kendall's tau to be computed.
The following are permitted:
\describe{
  \item{\code{MB}}{conditional Kendall's tau proposed in Martin and Betensky (2005) as \eqn{\hat{\tau}_c,}}
  \item{\code{IPW1}}{inverse probability weighted estimator proposed in Austin and Betensky (2014) as \eqn{\hat{\tau}_{c2},}}
  \item{\code{IPW2}}{restricted inverse probability weighted estimator proposed in Austin and Betensky (2014) as \eqn{\hat{\tau}_{c3}.}}
}}

\item{weights}{an optional vector of sampling weights used when \code{method = IPW1} or \code{method = IPW2}.
Inverse probability censored weighting (IPCW) is the default.}

\item{a}{a numeric transformation parameter. The default value is 0, which applies no transformation.
This parameter must be greater than \code{-1}.
See \code{?tranSurvfit} for the transformation model structure.}

\item{trans}{a character string specifying the transformation structure. The following are permitted:
\describe{
  \item{linear}{linear transformation structure,}
  \item{log}{log-linear transformation structure,}
  \item{exp}{exponential transformation structure.}
}}

\item{...}{for future methods.}
}
\value{
A numeric value representing the unconditional Kendall's tau.
}
\description{
Computes the conditional Kendall's tau and inference
}
\details{
This function performs statistical test for quasi-independence between truncation time and failure time.
The hypothesis test is based on the conditional Kendall's tau of Martin and Betensky (2005) and
the two versions of the inverse probability weighted Kendall's tau of Austin and Betensky (2014).

The output contains the following components:
\describe{
  \item{PE}{consistent point estimate of the conditional Kendall's tau.}
  \item{SE}{asymptotic standard error of the conditional Kendall's tau estimator.}
  \item{STAT}{the value of the normal test statistic.}
  \item{p.value}{the (Wald) p-value of the test.}
  \item{trans}{the transformation model (if applied).}
  \item{a}{the estimated transformation parameter.}
}
}
\examples{
data(channing, package = "boot")
chan <- subset(channing, sex == "Male" & entry < exit)
attach(chan)
cKendall(entry, exit, cens)
cKendall(entry, exit, cens, method = "IPW1")
cKendall(entry, exit, cens, method = "IPW2")
detach(chan)
}
\references{
Martin E. and Betensky R. A. (2005), Testing quasi-independence of failure and truncation times via conditional Kendall's tau,
\emph{Journal of the American Statistical Association}, \bold{100} (470): 484-492.

Austin, M. D. and Betensky R. A. (2014), Eliminating bias due to censoring in Kendall's tau estimators for quasi-independence of truncation and failure,
\emph{Computational Statistics & Data Analysis}, \bold{73}: 16-26.
}
\seealso{
\code{\link{trSurvfit}}
}
