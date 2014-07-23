## Suite of functions for testing gene expression amplitude


#' Test for differences in expression for one gene with a t-test
#'
#' @param log2cpm vector of log2-count-per-million values, the measure of expression
#' @param design design matrix for the linear model fit
#' @param test.variable the variable (or coefficient) for which testing is done
#' @return dataframe with logFC, SE, N (number of cells expressed), t-statistic and p-value
#' @seealso \code{\link{nchar}} which this function wraps
#' @export
#' @examples
#' x <- rnorm(10)
#' design <- model.matrix(~factor(rep(c(1,2), each = 5)) + factor(rep(c(1,2), 5)))
#' testAmp.ttest.onegene(x, design)
testAmp.ttest.onegene <- function(log2cpm, design, test.variable = ncol(design)) {
    fit <- lm(log2cpm ~ 0 + design)
    summ <- summary(fit)
    logFC <- SE <- t.value <- p.value <- NA
    if( test.variable <= nrow(summ$coef) ) {
        logFC <- summ$coef[test.variable, 1]
        SE <- summ$coef[test.variable, 2]
        t.value <- summ$coef[test.variable, 3]
        p.value <- summ$coef[test.variable, 4]
    }
    N <- length(log2cpm)
    data.frame(logFC, SE, N, t.value, p.value)
}
