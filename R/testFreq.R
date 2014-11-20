## Suite of functions for testing gene expression Frequency

#' Test for differences in expression frequency using logistic regression
#'
#' @param data_object a DGEList object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param design a design matrix to be used for the GLM fit.
#' @param coef the coefficient(s)
#' @param p_adj_method method to pass to \code{p.adjust()} to adjust
#' p-values for multiple testing.
#' @details Test for differences in expression frequency using
#' logistic regression and a likelihood ratio test. Returns a DGELRT
#' object with the differential frequency test results.
#' @export
#'
testFreq <- function(data_object, design, coef = ncol(design), p_adj_method = "BH", ntest = nrow(data_object)) {
    ## Define null design matrix
    design_null <- design[, -coef, drop = FALSE]
    ## Define matrix of binary "is expressed" values
    isexpr <- data_object$isexpr
    rownames(isexpr) <- rownames(data_object)
    ## Fit full model
    fit_full <- binomGLMFit(isexpr, design, ntest)
    ## Fit null model
    fit_null <- binomGLMFit(isexpr, design_null, ntest)
    ## Test for differential frequency
    LR <- fit_null$deviance - fit_full$deviance
    df_test <- fit_null$df_residual - fit_full$df_residual
    LRT_pvalue <- pchisq(LR, df = df_test, lower.tail = FALSE, log.p = FALSE)
    ## Compile results for output
    rn <- rownames(isexpr)
    if(is.null(rn))
        rn <- 1:nrow(isexpr)
    else
        rn <- make.unique(rn)
    tab <- data.frame(logFC = fit_full$coefficients[, coef], OverallFreq = rowMeans(isexpr),
                      LR = LR, PValue = LRT_pvalue, row.names = rn)
    fit_full$table <- tab
    fit_full$comparison <- colnames(design)[coef]
    fit_full$df_test <- df_test
    new("DGELRT", unclass(fit_full))
}


#' Test for differences in expression frequency using logistic regression
#'
#' @param y a matrix containing binary data on whether a
#' gene is expressed or not in each cell, with rows corresponding to
#' genes and columns to cells.
#' @param design a design matrix to be used for the GLM fit.
#' @param ntest the number of tests to carry out (default: test all genes).
#' @details Fit a binomial GLM for expression frequency to each gene.
#' @export
#'
binomGLMFit <- function(y, design, ntest = nrow(y)) {
    ## Set up fit object
    navec <- rep(NA, nrow(y))
    names(navec) <- rownames(y)
    namat <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
    rownames(namat) <- rownames(y)
    colnames(namat) <- colnames(y)
    fit <- list()
    fit$coefficients <- matrix(NA, nrow = nrow(y), ncol = ncol(design))
    colnames(fit$coefficients) <- colnames(design)
    rownames(fit$coefficients) <- rownames(y)
    fit$fitted_values <- fit$weights <- namat
    fit$isexpr <- y
    fit$rank <- fit$deviance <- fit$aic <- fit$null_deviance <- fit$iter <-
                fit$df_residual <- fit$df_null <- fit$converged <- fit$boundary <- navec
    ## Loop through genes getting binomial GLM fit
    for(i in 1:ntest) {
        y_gene <- y[i, ]
        fit_gene <- glm.fit(design, y_gene, family = binomial())
        fit$coefficients[i, ] <- fit_gene$coefficients
        fit$fitted_values[i, ] <- fit_gene$fitted.values
        fit$weights[i, ] <- fit_gene$weights
        fit$rank[i] <- fit_gene$rank
        fit$deviance[i] <- fit_gene$deviance
        fit$aic[i] <- fit_gene$aic
        fit$null_deviance[i] <- fit_gene$null.deviance
        fit$iter[i] <- fit_gene$iter
        fit$df_residual[i] <- fit_gene$df.residual
        fit$df_null[i] <- fit_gene$df.null
        fit$converged[i] <- fit_gene$converged
        fit$boundary[i] <- fit_gene$boundary
    }
    new("DGEGLM", fit)
}
