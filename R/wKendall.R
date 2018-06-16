#' Weighted conditional Kendall's tau
#'
#' This is a temporary function, could be merged into the main function \code{condKendall}.
#' Give conditional Kendall's tau with the ability to take on perturbation weights.
#'
#' @param trun left truncation time satisfying \code{trun} <= \code{obs}.
#' @param obs observed failure time, must be the same length as \code{trun}, might be right-censored.
#' @param weights perturbation weights.
#' 
#' @export
wKendall <- function(trun, obs, delta = NULL, weights = NULL) {
    n <- length(obs)
    if (is.null(delta)) delta <- rep(1, n)
    if (is.null(weights)) weights <- rep(1, n)
    .C("wKendallC", as.double(trun), as.double(obs), as.integer(n),
       as.double(delta), as.double(weights),
       out = double(1), PACKAGE = "tranSurv")$out
}
