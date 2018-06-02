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

