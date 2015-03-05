## Suite of functions for testing gene expression amplitude


#' Test for differences in expression amplitude
#'
#' @param object SCESet class object containing data
#' @param design design matrix describing experimental design
#' @param coef coefficient(s) of the design matrix for which testing is to be 
#' done
#' @param expr_values which expression values to use?
#' @param method method to be used to test for differences in expression 
#' amplitude (default is a t-test)
#' @param p_adj_method method to be used to adjust p-values for multiple testing
#' (defaults to Benjamini-Hochberg FDR)
#' @param ... arguments passed to \code{testAmp.ttest}
#' @details A user-level wrapper function for different approaches to
#' testing for differences in expression amplitude. Methods
#' implemented are: t-test, rank-production test, ...
#' @export
#'
testAmp <- function(object, design, coef = ncol(design), expr_values = "cpm",
                    method = "ttest", p_adj_method = "BH", ...) {
    method <- match.arg(method, c("ttest","rankprod"))
    if( method == "ttest" )
        results <- testAmp.ttest(object, expr_values, ...)
    else
        results <- "Sorry, other methods not implemented yet"
    results
}

#' Test for differences in expression amplitude using a t-test
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param expr_values character string giving the name of the element
#' of the SCESet object yielding the expression values to be used.
#' @param design design matrix describing experimental design
#' @param coef coefficient(s) of the design matrix for which testing is to be 
#' done
#' @param p_adj_method method to be used to adjust p-values for multiple testing
#' (defaults to Benjamini-Hochberg FDR)
#' @export
#'
testAmp.ttest <- function(object, design, coef = ncol(design), 
                          expr_values = "log2cpm.pc0.25", p_adj_method = "BH") {
    ## Define expression values
    if( is.null(object[[expr_values]]) )
        stop("Expression values not found")
    else
        expr <- object[[expr_values]]
    ## Define experiment and treatment vectors
    if( is.null(object$samples[[experiment]]) )
        stop("Experiment not found")
    else
        experiment <- object$samples[[experiment]]
    if( is.null(object$samples[[treatment]]) )
        stop("Treatment not found")
    else
        treatment <- object$samples[[treatment]]
    ## Check that we have informpation on whether genes in cells are expressed
    if( is.null(object$isexpr_cpm) )
        stop("Information on whether genes in cells are expressed not found. 
             Looked for object$isexpr_cpm.")
    ## Set up results object
    navec <- rep(NA, nrow(object))
    results <- data.frame(Ave.CPM = navec, logFC = navec, SE = navec, N = navec,
                          t.value = navec, p.value = navec)
    rownames(results) <- rownames(object)
    ## Conduct t-test for each gene
    for(i in 1:nrow(results)) {
        ## Code for one gene to fit lm
        dofit <- TRUE
        isexpr <- object$isexpr_cpm[i,]
        experiment.gene <- experiment[isexpr]
        treatment.gene <- treatment[isexpr]
        if( (length(unique(experiment.gene)) > 1) & 
                (length(unique(treatment.gene)) > 1) ) {
            design <- model.matrix(~ experiment.gene + treatment.gene)
            colnames(design) <- c("(Intercept)", "Expt", "Treatment")
        }
        else if( (length(unique(experiment.gene)) > 1) & 
                     (length(unique(treatment.gene)) <= 1) ) {
            warning(paste0("Gene ", rownames(results)[i], "expressed in only one
                           treatment. Cannot do test on treatment."))
            results$N[i] <- sum(isexpr)
            dofit <- FALSE
        } else if( length(unique(treatment.gene)) > 1 ) {
            warning(paste0("Gene ", rownames(results)[i], " expressed in only 
                           one experiment (experiment ", unique(experiment), "). 
                           Results should be treated with caution."))
            design <- model.matrix(~ treatment.gene)
            colnames(design) <- c("(Intercept)", "Metformin")
        } else {
            warning("Gene expressed in only one experiment and one treatment. 
                    No DE results possible.")
            results$N[i] <- sum(isexpr)
            dofit <- FALSE
        }
        if(dofit) {
            results$Ave.CPM[i] <- mean(object$log2cpm.pc0.25[i, isexpr])
            fit <- testAmp.ttest.onegene(object$log2cpm.pc0.25[i, isexpr], 
                                         design)
            results[i, -1] <- fit
        }
    }
    ## Adjust p-values for multiple testing
    p.adj.method <- match.arg(p.adj.method,
                              c("holm", "hochberg", "hommel", "bonferroni", 
                                "BH", "BY", "fdr", "none"))
    results$FDR <- p.adjust(results$p.value, method = p.adj.method)
    ## Output results
    results
}


#' Test for differences in expression for one gene with a t-test
#'
#' @param log2cpm vector of log2-count-per-million values, the measure of 
#' expression
#' @param design design matrix for the linear model fit
#' @param test_variable the variable (or coefficient) for which testing is done
#' @return dataframe with logFC, SE, N (number of cells expressed), t-statistic 
#' and p-value
#' @seealso \code{\link{nchar}} which this function wraps
#' @export
#' @examples
#' x <- rnorm(10)
#' design <- model.matrix(~factor(rep(c(1,2), each = 5)) + factor(rep(c(1,2), 5)))
#' testAmp.ttest.onegene(x, design)
testAmp.ttest.onegene <- function(log2cpm, design, test_variable = ncol(design)) {
    fit <- lm(log2cpm ~ 0 + design)
    summ <- summary(fit)
    logFC <- SE <- t.value <- p.value <- NA
    if( test_variable <= nrow(summ$coef) ) {
        logFC <- summ$coef[test_variable, 1]
        SE <- summ$coef[test_variable, 2]
        t.value <- summ$coef[test_variable, 3]
        p.value <- summ$coef[test_variable, 4]
    }
    N <- length(log2cpm)
    data.frame(logFC, SE, N, t.value, p.value)
}

