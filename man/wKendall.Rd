% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wKendall.R
\name{wKendall}
\alias{wKendall}
\title{Weighted conditional Kendall's tau}
\usage{
wKendall(trun, obs, delta = NULL, weights = NULL)
}
\arguments{
\item{trun}{left truncation time satisfying \code{trun} <= \code{obs}.}

\item{obs}{observed failure time, must be the same length as \code{trun}, might be right-censored.}

\item{delta}{an optional 0-1 vector of censoring indicator (0 = censored, 1 = event) for \code{obs}.
If this vector is not specified, \code{cKendall} assumes no censoring and all observed failure time
denote events.}

\item{weights}{an optional perturbation weights.}
}
\value{
A numeric value representing the weighted conditional Kendall's tau.
}
\description{
This is function computes the perturbed version of the conditional Kendall's tau.
}
\examples{
data(channing, package = "boot")
chan <- subset(channing, sex == "Male" & entry < exit)
## When weights is not specified, this function reduces to condKendall()
with(chan, wKendall(entry, exit, cens))
mean(replicate(1000, with(chan, wKendall(entry, exit, cens, rexp(nrow(chan))))))
}
