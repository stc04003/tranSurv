#' @export
print.condKendall <- function(x, ...) {
    cat("\n Test for quasi-independence with conditional Kendall's tau\n")
    cat("\n Call: ")
    print(x$Call)
    if (x$a != 0)
        cat(paste("\nTransformation is applied with parameter a =", round(x$a, 4)))
    cat(paste("\n", "Kendall's tau =", round(x$PE, 4), ", SE =", round(x$SE, 4),
              ", Z =", round(x$STAT, 4), ", p-value = ", round(x$p.value, 4), "\n\n"))
}

#' @export
print.pmcc <- function(x, ...) {
    cat("\n Test for quasi-independence with conditional correlation coefficient\n")
    cat("\n Call: ")
    print(x$Call)
    if (x$a != 0)
        cat(paste("\nTransformation is applied with parameter a =", round(x$a, 4)))
    cat(paste("\n", "Correlation coefficient =", round(x$PE, 4), ", SE =", round(x$SE, 4),
              ", Z =", round(x$STAT, 4), ", p-value = ", round(x$p.value, 4), "\n\n"))
}

#' @export
print.trSurvfit <- function(x, ...) {
    cat("\n Fitting structural transformation model \n")
    cat("\n Call: ")
    print(x$Call)
    cat(paste("\n", "Conditional Kendall's tau =",
              round(x$iniKendall, 4), ", p-value =", round(x$iniP, 4)))
    cat(paste("\n", "Restricted inverse probability weighted Kendall's tau =",
              round(x$iniKendall.ipw, 4), ", p-value =", round(x$iniP.ipw, 4)))
    cat(paste("\n Transformation parameter by minimizing absolute value of Kendall's tau:",
              round(x$byTau$par, 4)))
    cat(paste("\n Transformation parameter by maximizing p-value of the test:",
              round(x$byP$par, 4), "\n\n"))
}

#' @export
#' @importFrom stats model.matrix printCoefmat sd
print.trReg <- function(x, ...) {
    cat("\n Call:")
    print(x$Call)
    cat("\n Sample size =", nrow(x$.data))
    cat("\n Number of events = ", sum(x$.data$status))
    cat("\n Transformation parameter is", x$a, "\n\n")
    tab <- cbind(coef = round(x$PE[,1], 3),
                 "se(coef)" = round(x$SE, 3),
                 z = round(x$PE[,1] / x$SE, 3),
                 "Pr(>|z|)" = round(2 * pnorm(-abs(x$PE[,1] / x$SE)), 3))
    rownames(tab) <- x$varNames
    printCoefmat(as.data.frame(tab), P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
}

#' @export
print.trgof <- function(x, ...) {
    cat("\n Overall signficances based on left-truncated regression model: p-value =", round(x$pval))
    cat("\n Transformation model based on", x$Q, "segements.")
    vv <- NULL
    for (i in 1:(length(x$breaks) - 1)) vv[i] <- paste("(", x$breaks[i], ", ", x$breaks[i + 1], "]", sep = "")
    cat("\n The segments are", paste(vv, collapse = " "), "\n\n")
    if (!is.null(x$fit.all)) {
        printCoefmat(round(coef(summary(x$fit.all)), 4), P.values = TRUE, has.Pvalue = TRUE)
        cat("\n")
    }
}
