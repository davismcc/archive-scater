## all code related to the RankProd methods
## Davis McCarthy
## 18 November 2015

##----------------------------------------------------------------------------##
# 
# # Compare and aggregate results from different ranking methods using exact rank-product tests
# # Davis McCarthy, July 2014
# 
# #' Aggregate and compare results from different ranking methods using exact rank-product tests
# #'
# #' @param results a matrix or data.frame object where columns provide
# #' results values (e.g. p-values) to be used for ranking features
# #' (e.g. genes) and ordering of results.
# #' @param nmax_genes integer giving the maximum number of features
# #' (e.g. genes) for which to compute exact rank-product p-values
# #' @param max_rankprod_exact integer giving the maximum rank-product for
# #' which to compute exact p-values. Exact p-values are very slow to
# #' compute for large rank-products, and converge to approximate
# #' p-values.
# #' @return data frame with mean rank, rank product, number of ways of
# #' obtaining given rank product, approximate p-value and exact p-value
# #' @import matrixStats
# #' @export
# aggregateResults <- function(results, nmax_genes = 500, max_rankprod_exact = 100000) {
#     ranks <- matrix(NA, nrow = nrow(results), ncol = ncol(results), dimnames = list(rownames(results), colnames(results)))
#     ## Define ranking for each feature, for each method
#     for(i in seq_len(ncol(results))) {
#         ranks[, i] <- rank(results[, i], ties.method = "average")
#     }
#     rank_products <- matrixStats::rowProds(ranks)
#     names(rank_products) <- rownames(ranks)
#     ranks <- ranks[order(rank_products), ]
#     pvals <- calcRankProdPValsBounds(ranks, n = length(ranks), k=ncol(results))
#     aggregated_results <- data.frame(Feature=rownames(results), 
#                                      Rank_Prod=rank_products, P_Value=pvals, 
#                                      Rank=ranks)
#     aggregated_results[1:nmax_genes, ]
# }
# 
# ##----------------------------------------------------------------------------##
# 
# #' Right-tailed gamma approximation to rank-product p-value
# #'
# #' @param r integer rank-product
# #' @param k integer number of experiments
# #' @param n integer number of genes
# #' @return float aproximate p-value
# #' @export
# #' @examples
# #' calcRankProdPvalsApprox(9720, 5, 500)
# calcRankProdPvalsApprox <- function(r, k, n) {
#     1 - pgamma(-log(r/(n+1)^k), k, scale = 1)
# }
# 
# 
# #' Calculate bounds on p-values for rank-product statistics
# #'
# #'
# #' This function computes bounds on the p-value for rank products.
# #'
# #'
# #' @param rho a vector of integers corresponding to the rank products
# #' for which one wishes to compute the p-value.
# #' @param n the number of features (e.g. genes).
# #' @param k the number of replicates.
# #' @param delta a character string indicating whether an upper bound
# #' ('upper'), lower bound ('lower'), or the default,geometric
# #' approximation (i.e. the geometric mean of the upper and lower
# #' bounds) ('geometric') should be computed.
# #'
# #' @return a vector of p-values, one for each rank product.
# #'
# #' @details The exact p-value is guaranteed to be in between the lower
# #' and the upper bound. The geometric mean of the two bounds can be
# #' used as an approximation. Each bound is a piecewise continuous
# #' function of the rank product. The different pieces each have an
# #' analytic form, the parameters of which can be computed recursively.
# #'
# #' This implementation closely follows the description in Heskes,
# #' Eisinga, Breitling: "A fast algorithm for determining bounds and
# #' accurate approximate p-values of the rank product statistic for
# #' replicate experiments", further referred to as HEB. More
# #' specifically, this R function corresponds to the recursive variant,
# #' sketched as pseudocode in the additional material of HEB.
# #'
# #' We thank HEB for making their R code available. Their original code
# #' has been lightly modified for use in \pkg{scater}.
# #'
# #' @export
# #' @examples
# #' calcRankProdPValsBounds(c(10, 100, 1000), 2000, 10, "upper")
# calcRankProdPValsBounds <- function(rho, n, k, delta = "geometric") {
#     ## INPUT HANDLING
#     if(any(rho > n^k) || any(rho < 1))
#         stop('rho out of bounds')
#     if(is.numeric(delta) == FALSE) {
#         if(delta == 'geometric') {
#             temp1 <- calcRankProdPValsBounds(rho, n, k, 'upper')
#             temp2 <- calcRankProdPValsBounds(rho, n, k, 'lower')
#             pvalue <- sqrt(temp1*temp2)   # geometric mean of upper and lower bound
#             return(pvalue)
#         }
#         else {
#             delta <- switch(delta,
#                             upper = 1,        # for computing upper bound
#                             lower = 0)        # for computing lower bound
#         }
#     }
#     ## COMPUTE INTERVALS THAT CONTAIN THE RANK PRODUCTS
#     logn <- log(n)
#     allj <- ceiling(-(log(rho) / logn) + k)   # index specifying the interval that contains rho 
#     minj <- min(allj)                     # lowest interval index
#     maxj <- max(allj)                     # highest interval index
#     ## INITIALIZE PARAMETERS
#     param <- matrix(list(), nrow = k+1, ncol = maxj+1)
#     for(i in 1:(k + 1)){
#         for(j in 1:(maxj + 1)){
#             param[[i, j]] <- list(a = c(), b = c(), c = c(), d = c(), e = c())
#         }
#     }
#     ## param is a matrix of lists; each element of param is a list with values for the parameters
#     ## a through e, which correspond to the parameters alpha through epsilon in HEB;
#     ## specifially, param[[i+1,j+1]]$a corresponds to alpha_{i,j} in HEB, etc, where the offset
#     ## of 1 is introduced to be able to represent, for example, alpha_{0,0};
#     ## a, b, and c can be vectors (with possibly different lengths for different i and j),
#     ## d and e are scalars
#     ## COMPUTE PARAMETERS
#     for(j in minj:maxj) {
#         param <- updateParam(param, n, k, j, delta)
#     }
#     ## call to the function updateParam which recursively computes all parameters that are needed
#     ## to calculate the p-value for a rank product rho that lies in the interval with index j
#     ## COMPUTE RANK PRODUCTS GIVEN PARAMETERS
#     k1 <- 1+k
#     G <- rep(0, length(rho))   # G is a vector of the same length as rho,
#     # for each rho bounding the number of rank products 
#     for(j in minj:maxj) {
#         j1 <- 1+j
#         iii <- which(allj == j)         # indices of all rank products that fall in interval j:
#         # bounds for these rank products can be computed with
#         # the same set of parameters                                    
#         thisrho <- rho[iii]
#         thisparam <- param[[k1, j1]]
#         thisG <- thisparam$e
#         if(j != 0) {
#             nrho <- length(thisrho)
#             nterms <- length(thisparam$a)
#             thisG <- thisG + thisparam$d * thisrho
#             d1 <- matrix(thisparam$c) %*% thisrho
#             d2 <- matrix(rep(log(thisrho), nterms), nrow = nterms, byrow = TRUE) -
#                 t(matrix(rep(logn * (k - j + thisparam$b), nrho), nrow = nrho, byrow = TRUE))
#             d3 <- t(matrix(rep(thisparam$a, nrho), nrow = nrho, byrow = TRUE)) 
#             thisG <- thisG + colSums(d1 * (d2^d3))
#         }
#         # the 10 lines above implement equation (8) in HEB
#         G[iii] <- thisG
#     }
#     pvalue <- G/n^k
#     return(pvalue)
# }
# 
# #' Update parameters for computation of rank-product p-values
# #'
# #' This subroutine updates the current set of parameters to make sure
# #' that the parameters corresponding to k replicates and the j'th
# #' interval are included.
# #'
# #' @param param a matrix of lists, where each element of param is a
# #' list with values for the parameters a through e; these parameters
# #' specify the functional form of the bound; a, b, and c are all
# #' vectors of unknown length, d and e are scalars.
# #' @param n the number of features (e.g. genes).
# #' @param k the number of replicates for which we need to compute the
# #' corresponding parameters.
# #' @param j the index of the interval for which we need to compute the
# #' corresponding parameters.
# #' @param delta 0 for the lower bound and 1 for the upper bound.
# #'
# #' @return A possibly updated set of parameters, at least including
# #' those corresponding to (k,j).
# #'
# #' @details This subroutine makes sure that the parameters
# #' corresponding to k replicates and a rank product within the j'th
# #' interval are already present. If they already are (because
# #' calculated before), it does not compute anything. Otherwise, it
# #' recursively computes all parameters that are needed to arrive at
# #' the parameters for (k,j).
# #'
# #' This implementation closely follows HEB, in particular equations (9) through (11).
# #' 
# updateParam <- function(param, n, k, j, delta) {
#     k1 <- 1+k
#     j1 <- 1+j
#     if(length(param[[k1,j1]]$e) == 0) {  # apparently empty, so needs to be calculated
#         if(j == 0) {   # initializing G_{k0}
#             param[[k1,j1]]$e <- n^k
#             param[[k1,j1]]$d <- 0
#             # the 2 lines above implement equation (11) in HEB
#         }
#         else {
#             k0 <- k1-1
#             j0 <- j1-1
#             param <- updateParam(param, n, k-1, j-1, delta)
#             # checking that the parameters for (k-1,j-1) that are needed to compute the
#             # parameters for (k,j) are indeed available; if not, they are themselves computed
#             param00 = param[[k0, j0]]
#             newa0 = param00$a + 1
#             newb0 = param00$b
#             newc0 = param00$c / newa0
#             param11 = param00
#             # the 5 lines above predefine some parameters common to equations (9) and (10) in HEB
#             if(k == j){ # updates for G_{kk}
#                 param11$e <- (1 - delta) * (1 - param00$e)
#                 param11$d <- delta * param00$d + param00$e
#                 param11$a <- c(1, param00$a, newa0)
#                 param11$b <- c(0, param00$b, newb0)
#                 param11$c <- c(param00$d, delta * param00$c, newc0)
#                 # the 5 lines above implement equation (10) in HEB
#             }
#             else {  # updates for G_{kj}, j < k
#                 param <- updateParam(param, n, k-1, j, delta)
#                 # checking that the parameters for (k-1,j) that are needed to compute the
#                 # parameters for (k,j) are indeed available; if not, they are themselves computed
#                 param01 <- param[[k0, j1]]
#                 logn <- log(n)
#                 lognnkj <- (k - j) * logn
#                 newa1 <- param01$a + 1
#                 newa <- c(newa0, newa1)
#                 newb <- c(newb0, param01$b)
#                 newc <- c(newc0, -param01$c / newa1)
#                 param11$e <- n * param01$e + (delta - 1) * (param00$e - param01$e)
#                 lognminb <- c(-1 * param00$b * logn, (1 - param01$b) * logn)
#                 param11$d <- delta * param00$d + (1 - delta) * param01$d / n + 
#                     (param00$e - param01$e) / exp(lognnkj) - sum(newc * (lognminb^newa))
#                 param11$a <- c(1, 1, param00$a, param01$a, newa)
#                 param11$b <- c(0, 1, param00$b, param01$b, newb)
#                 param11$c <- c(param00$d, -param01$d, delta * param00$c, (1-delta) * param01$c / n, newc)
#                 # the 15 lines above implement equation (9) in HEB
#             }
#             param[[k1, j1]] <- makeUnique(param11)
#             # although not strictly necessary, the a, b and c vectors can possibly be shortened by
#             # restricting oneselves to unique combinations of a and b values
#         }
#     }
#     return(param)
# }
# 
# #' Update the parameters so that they are unique
# #'
# #' This subroutine updates the parameters for a specific number of replicates and interval
# #' such that it contains only unique combinations of the parameters a and b.
# #'
# #' @param param a single list with values for the parameters a through
# #' e; these parameters specify the functional form of the bound; a, b,
# #' and c are all vectors of unknown length, d and e are scalars.
# #' 
# #' @return A possibly updated and then more concise set of parameters
# #' containing only unique combinations of the parameters a and b.
# #'
# #' @details While updating the vectors a and b, one may end up with
# #' the exact same combinations of a and b. Given the functional form
# #' of the bound, the representation can then be made more concise by
# #' simply adding the corresponding elements of c.
# makeUnique <- function(param) {
#     ab <- t(rbind(param$a, param$b))
#     uniqueab <- unique(ab)
#     nunique <- dim(uniqueab)[1]
#     param$a <- t(uniqueab[, 1])
#     param$b <- t(uniqueab[, 2])
#     newc <- rep(0, nunique)
#     for(i in 1:nunique) {
#         iii <- intersect(which(ab[, 1] == uniqueab[i, 1]), which(ab[, 2] == uniqueab[i, 2]))
#         newc[i] <- sum(param$c[iii])  
#     }
#     param$c <- newc
#     return(param)
# }
# 
# ##----------------------------------------------------------------------------##
# 
# #' Get logFC values for replicates within one experimental block
# #'
# #' @param data_matrix numeric matrix with expression values
# #' (log-scale) for features (rows) and samples/cells (columns).
# #' @param group character vector or factor defining the two groups to
# #' be compared
# #' @return matrix of logFC values for all pairwise comparisons of
# #' samples in group1 and group2
# #' @export
# getLogFCOneBlock <- function(data_matrix, group) {
#     ## Check that group has the correct number of levels: only two groups currently supported
#     group <- as.factor(group)
#     if(nlevels(group) > 2)
#         stop("More than two groups to be compared. Can only compare between two groups, currently.")
#     ## Define indices for groups 1 and 2
#     idx1 <- which(group == levels(group)[1])
#     idx2 <- which(group == levels(group)[2])
#     ngroup1 <- length(idx1)
#     ngroup2 <- length(idx2)
#     ## Define number of genes, comparisons and FC matrix
#     ngenes <- nrow(data_matrix)
#     ncells <- ncol(data_matrix)
#     k <- ngroup1 * ngroup2
#     logFC <- matrix(NA, nrow = ngenes, ncol = k)
#     ## Get logFC for all pairwise comparisons of genes in group1 and group2
#     for(i in seq_len(ngroup1)) {
#         idx_this_insertion <- ((i - 1) * ngroup2 + 1):(i * ngroup2)
#         insert_this <- data_matrix[, idx1[i]] - data_matrix[, idx2]
#         logFC[, idx_this_insertion] <- insert_this
#     }
#     ## Return logFC matrix
#     logFC
# }
# 
# ##----------------------------------------------------------------------------##
# 
# #' Get ranks from logFCs for one block
# #'
# #' Get ranks both for up-regulation in group1 and up-regulation in
# #' group2
# #'
# #' @param logFC_matrix numeric matrix of logFCs computed as group1 -
# #' group2
# #' @param group character vector or factor defining the experimental
# #' groups being compared
# #' @param logged boolean: return raw ranks if \code{FALSE} or
# #' log(ranks) if \code{TRUE}
# #' @export
# getRanksOneBlock <- function(logFC_matrix, group, logged = FALSE) {
#     group <- as.factor(group)
#     label1 <- paste0(levels(group), collapse = "-")
#     label2 <- paste0(rev(levels(group)), collapse = "-")
#     ## We assume that logFCs have been computed such that group1 = levels(group)[1] and logFC is group1 - group2
#     ## Ranks in ascending order, so make logFC negative so that large positive logFCs end up at the top of the ranking
#     rank_matrix_up_grp1 <- apply(-logFC_matrix, 2, rank)
#     ## Rank by up-regulation in group2
#     rank_matrix_up_grp2 <- apply(logFC_matrix, 2, rank)
#     if(logged) {
#         rank_matrix_up_grp1 <- log(rank_matrix_up_grp1)
#         rank_matrix_up_grp2 <- log(rank_matrix_up_grp2)
#     }
#     ranks <- list(rank_matrix_up_grp1, rank_matrix_up_grp2)
#     names(ranks) <- c(label1, label2)
#     ranks
# }
# 
# ##----------------------------------------------------------------------------##
# 
# ## Scale expression values to have average amplitude of 1 in
# ## expressing cells
# 
# #' Scale expression values to have average amplitude of 1 in
# #' expressing cells
# #'
# #' @param expr_matrix a numeric matrix containing expression values
# #' counts.
# #' @param is_expr a logical matrix of same dimensions as
# #' \code{expr_matrix} indicating which cells are expressed for each
# #' gene.
# #' @details This function is intended to process RNA-Seq or ChIP-Seq
# #' data prior to modeling with RankProduct. Written by Davis McCarthy.
# #'
# #' @return A matrix object.
# #' @export
# #'
# scaleAmp <- function(expr_matrix, is_expr) {
#     ## Get average amplitude in expressing cells for each gene
#     expr_vals <- expr_matrix
#     expr_vals[!is_expr] <- NA
#     ave_expr <- rowMeans(expr_vals, na.rm = TRUE)
#     ## ave_expr <- rep(NA, nrow(expr_matrix))
#     ## for( i in seq_len(nrow(expr_matrix)) ) {
#     ##     expr_ind <- is_expr[i, ]
#     ##     expr_val <- expr_matrix[i, ]
#     ##     ave_expr[i] <- mean(expr_val[expr_ind], na.rm = TRUE)
#     ## }
#     ## Scale expression values for expressed observations
#     out <- expr_matrix / ave_expr
#     ## Leave expression values for non-expressed observations unscaled
#     ## out[!is_expr] <- expr_matrix[!is_expr]
#     out
# }
# 
# 
# ## Scale expression values to have unit variance in expressing cells
# 
# #' Scale expression values to have unit variance in expressing cells
# #'
# #' @param expr_matrix a numeric matrix containing expression values
# #' counts.
# #' @param is_expr a logical matrix of same dimensions as
# #' \code{expr_matrix} indicating which cells are expressed for each
# #' gene.
# #' @param vars optional vector of variance values to apply for the
# #' scaling of the expression values. If \code{NULL} (default), then
# #' variance is computed for each gene only from the cells with
# #' non-zero expression for the gene.
# #' @details This function is intended to process RNA-Seq or ChIP-Seq
# #' data prior to modeling with RankProduct. Written by Davis McCarthy.
# #'
# #' @return A matrix object.
# #' @export
# #'
# scaleVar <- function(expr_matrix, is_expr, vars = NULL) {
#     if( is.null(vars) ) {
#         ## Get variance of expression values in expressing cells for each gene
#         ## Use rowVars() from the matrixStats package
#         expr_vals <- expr_matrix
#         expr_vals[!is_expr] <- NA
#         var_expr <- matrixStats::rowVars(expr_vals, na.rm = TRUE)
#         ## If any variances are zero, convert to 1 to keep observations as is
#         var_expr[var_expr == 0] <- 1
#         ## If any variances are NA, similarly replace the NA with 1
#         var_expr[is.na(var_expr)] <- 1
#     }
#     else {
#         var_expr <- vars
#     }
#     ## Scale expression values for expressed observations
#     out <- expr_matrix / sqrt(var_expr)
#     ## Leave expression values for non-expressed observations unscaled
#     ## out[!is_expr] <- expr_matrix[!is_expr]
#     out
# }
# 
# ##----------------------------------------------------------------------------##
# 
# ## Suite of functions for testing gene expression amplitude
# 
# 
# #' Test for differences in expression amplitude
# #'
# #' @param object SCESet class object containing data
# #' @param design design matrix describing experimental design
# #' @param coef coefficient(s) of the design matrix for which testing is to be 
# #' done
# #' @param expr_values which expression values to use?
# #' @param method method to be used to test for differences in expression 
# #' amplitude (default is a t-test)
# #' @param p_adj_method method to be used to adjust p-values for multiple testing
# #' (defaults to Benjamini-Hochberg FDR)
# #' @param ... arguments passed to \code{testAmp.ttest}
# #' @details A user-level wrapper function for different approaches to
# #' testing for differences in expression amplitude. Methods
# #' implemented are: t-test, rank-production test, ...
# #' @export
# #'
# testAmp <- function(object, design, coef = ncol(design), expr_values = "cpm",
#                     method = "ttest", p_adj_method = "BH", ...) {
#     method <- match.arg(method, c("ttest","rankprod"))
#     if( method == "ttest" )
#         results <- testAmp.ttest(object, expr_values, ...)
#     else
#         results <- "Sorry, other methods not implemented yet"
#     results
# }
# 
# #' Test for differences in expression amplitude using a t-test
# #'
# #' @param object an SCESet object containing expression values and
# #' experimental information. Must have been appropriately prepared.
# #' @param expr_values character string giving the name of the element
# #' of the SCESet object yielding the expression values to be used.
# #' @param design design matrix describing experimental design
# #' @param coef coefficient(s) of the design matrix for which testing is to be 
# #' done
# #' @param p_adj_method method to be used to adjust p-values for multiple testing
# #' (defaults to Benjamini-Hochberg FDR)
# #' @export
# #'
# testAmp.ttest <- function(object, design, coef = ncol(design), 
#                           expr_values = "log2cpm.pc0.25", p_adj_method = "BH") {
#     ## Define expression values
#     if( is.null(object[[expr_values]]) )
#         stop("Expression values not found")
#     else
#         expr <- object[[expr_values]]
#     ## Define experiment and treatment vectors
#     if( is.null(object$samples[[experiment]]) )
#         stop("Experiment not found")
#     else
#         experiment <- object$samples[[experiment]]
#     if( is.null(object$samples[[treatment]]) )
#         stop("Treatment not found")
#     else
#         treatment <- object$samples[[treatment]]
#     ## Check that we have informpation on whether genes in cells are expressed
#     if( is.null(object$isexpr_cpm) )
#         stop("Information on whether genes in cells are expressed not found. 
#              Looked for object$isexpr_cpm.")
#     ## Set up results object
#     navec <- rep(NA, nrow(object))
#     results <- data.frame(Ave.CPM = navec, logFC = navec, SE = navec, N = navec,
#                           t.value = navec, p.value = navec)
#     rownames(results) <- rownames(object)
#     ## Conduct t-test for each gene
#     for(i in 1:nrow(results)) {
#         ## Code for one gene to fit lm
#         dofit <- TRUE
#         isexpr <- object$isexpr_cpm[i,]
#         experiment.gene <- experiment[isexpr]
#         treatment.gene <- treatment[isexpr]
#         if( (length(unique(experiment.gene)) > 1) & 
#             (length(unique(treatment.gene)) > 1) ) {
#             design <- model.matrix(~ experiment.gene + treatment.gene)
#             colnames(design) <- c("(Intercept)", "Expt", "Treatment")
#         }
#         else if( (length(unique(experiment.gene)) > 1) & 
#                  (length(unique(treatment.gene)) <= 1) ) {
#             warning(paste0("Gene ", rownames(results)[i], "expressed in only one
#                            treatment. Cannot do test on treatment."))
#             results$N[i] <- sum(isexpr)
#             dofit <- FALSE
#         } else if( length(unique(treatment.gene)) > 1 ) {
#             warning(paste0("Gene ", rownames(results)[i], " expressed in only 
#                            one experiment (experiment ", unique(experiment), "). 
#                            Results should be treated with caution."))
#             design <- model.matrix(~ treatment.gene)
#             colnames(design) <- c("(Intercept)", "Metformin")
#         } else {
#             warning("Gene expressed in only one experiment and one treatment. 
#                     No DE results possible.")
#             results$N[i] <- sum(isexpr)
#             dofit <- FALSE
#         }
#         if(dofit) {
#             results$Ave.CPM[i] <- mean(object$log2cpm.pc0.25[i, isexpr])
#             fit <- testAmp.ttest.onegene(object$log2cpm.pc0.25[i, isexpr], 
#                                          design)
#             results[i, -1] <- fit
#         }
#     }
#     ## Adjust p-values for multiple testing
#     p.adj.method <- match.arg(p.adj.method,
#                               c("holm", "hochberg", "hommel", "bonferroni", 
#                                 "BH", "BY", "fdr", "none"))
#     results$FDR <- p.adjust(results$p.value, method = p.adj.method)
#     ## Output results
#     results
# }
# 
# 
# #' Test for differences in expression for one gene with a t-test
# #'
# #' @param log2cpm vector of log2-count-per-million values, the measure of 
# #' expression
# #' @param design design matrix for the linear model fit
# #' @param test_variable the variable (or coefficient) for which testing is done
# #' @return dataframe with logFC, SE, N (number of cells expressed), t-statistic 
# #' and p-value
# #' @seealso \code{\link{nchar}} which this function wraps
# #' @export
# #' @examples
# #' x <- rnorm(10)
# #' design <- model.matrix(~factor(rep(c(1,2), each = 5)) + factor(rep(c(1,2), 5)))
# #' testAmp.ttest.onegene(x, design)
# testAmp.ttest.onegene <- function(log2cpm, design, test_variable = ncol(design)) {
#     fit <- lm(log2cpm ~ 0 + design)
#     summ <- summary(fit)
#     logFC <- SE <- t.value <- p.value <- NA
#     if( test_variable <= nrow(summ$coef) ) {
#         logFC <- summ$coef[test_variable, 1]
#         SE <- summ$coef[test_variable, 2]
#         t.value <- summ$coef[test_variable, 3]
#         p.value <- summ$coef[test_variable, 4]
#     }
#     N <- length(log2cpm)
#     data.frame(logFC, SE, N, t.value, p.value)
# }
# 
# ##----------------------------------------------------------------------------##
# 
# ## Suite of functions for testing gene expression Frequency
# 
# #' Test for differences in expression frequency using logistic regression
# #'
# #' @param data_object a DGEList object containing expression values and
# #' experimental information. Must have been appropriately prepared.
# #' @param design a design matrix to be used for the GLM fit.
# #' @param coef the coefficient(s)
# #' @param p_adj_method method to pass to \code{p.adjust()} to adjust
# #' p-values for multiple testing.
# #' @param ntest number of tests to carry out (defaults to number of genes)
# #' @details Test for differences in expression frequency using
# #' logistic regression and a likelihood ratio test. Returns a DGELRT
# #' object with the differential frequency test results.
# #' @export
# #'
# testFreq <- function(data_object, design, coef = ncol(design), 
#                      p_adj_method = "BH", ntest = nrow(data_object)) {
#     ## Define null design matrix
#     design_null <- design[, -coef, drop = FALSE]
#     ## Define matrix of binary "is expressed" values
#     isexpr <- data_object$isexpr
#     rownames(isexpr) <- rownames(data_object)
#     ## Fit full model
#     fit_full <- binomGLMFit(isexpr, design, ntest)
#     ## Fit null model
#     fit_null <- binomGLMFit(isexpr, design_null, ntest)
#     ## Test for differential frequency
#     LR <- fit_null$deviance - fit_full$deviance
#     df_test <- fit_null$df_residual - fit_full$df_residual
#     LRT_pvalue <- pchisq(LR, df = df_test, lower.tail = FALSE, log.p = FALSE)
#     ## Compile results for output
#     rn <- rownames(isexpr)
#     if(is.null(rn))
#         rn <- 1:nrow(isexpr)
#     else
#         rn <- make.unique(rn)
#     tab <- data.frame(logFC = fit_full$coefficients[, coef], OverallFreq = rowMeans(isexpr),
#                       LR = LR, PValue = LRT_pvalue, row.names = rn)
#     fit_full$table <- tab
#     fit_full$comparison <- colnames(design)[coef]
#     fit_full$df_test <- df_test
#     new("DGELRT", unclass(fit_full))
# }
# 
# 
# #' Test for differences in expression frequency using logistic regression
# #'
# #' @param y a matrix containing binary data on whether a
# #' gene is expressed or not in each cell, with rows corresponding to
# #' genes and columns to cells.
# #' @param design a design matrix to be used for the GLM fit.
# #' @param ntest the number of tests to carry out (default: test all genes).
# #' @details Fit a binomial GLM for expression frequency to each gene.
# #' @export
# #'
# binomGLMFit <- function(y, design, ntest = nrow(y)) {
#     ## Set up fit object
#     navec <- rep(NA, nrow(y))
#     names(navec) <- rownames(y)
#     namat <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
#     rownames(namat) <- rownames(y)
#     colnames(namat) <- colnames(y)
#     fit <- list()
#     fit$coefficients <- matrix(NA, nrow = nrow(y), ncol = ncol(design))
#     colnames(fit$coefficients) <- colnames(design)
#     rownames(fit$coefficients) <- rownames(y)
#     fit$fitted_values <- fit$weights <- namat
#     fit$isexpr <- y
#     fit$rank <- fit$deviance <- fit$aic <- fit$null_deviance <- fit$iter <-
#         fit$df_residual <- fit$df_null <- fit$converged <- fit$boundary <- navec
#     ## Loop through genes getting binomial GLM fit
#     for(i in 1:ntest) {
#         y_gene <- y[i, ]
#         fit_gene <- glm.fit(design, y_gene, family = binomial())
#         fit$coefficients[i, ] <- fit_gene$coefficients
#         fit$fitted_values[i, ] <- fit_gene$fitted.values
#         fit$weights[i, ] <- fit_gene$weights
#         fit$rank[i] <- fit_gene$rank
#         fit$deviance[i] <- fit_gene$deviance
#         fit$aic[i] <- fit_gene$aic
#         fit$null_deviance[i] <- fit_gene$null.deviance
#         fit$iter[i] <- fit_gene$iter
#         fit$df_residual[i] <- fit_gene$df.residual
#         fit$df_null[i] <- fit_gene$df.null
#         fit$converged[i] <- fit_gene$converged
#         fit$boundary[i] <- fit_gene$boundary
#     }
#     new("DGEGLM", fit)
# }
# 
# ##----------------------------------------------------------------------------##
# 
# # Wrapper functions for the computation of p-values for rank-product tests
# # Davis McCarthy, July 2014
# 
# #' Conduct rank-product test
# #'
# #' @param dfilt DGEList object containing data and so on that has been
# #' filtered and expressed cells determined
# #' @param nmax.genes integer giving the maximum number of genes for
# #' which to compute exact rank-product p-values
# #' @param max.rprod.exact integer giving the maximum rank-product for
# #' which to compute exact p-values. Exact p-values are very slow to
# #' compute for large rank-products, and converge to approximate
# #' p-values.
# #' @return data frame with mean rank, rank product, number of ways of
# #' obtaining given rank product, approximate p-value and exact p-value
# #' 
# testRankProduct <- function(dfilt, nmax.genes = 500, max.rprod.exact = 100000) {
#     
#     #### NEED TO ADJUST TO GENERIC FOR SCESet Class Object #####
#     
#     #     navec <- rep(NA, nrow(dfilt))
#     #     rp.results <- data.frame(Ave.CPM = navec, Ave.logFC = navec, N = navec)
#     #     rownames(rp.results) <- rownames(dfilt)
#     #     rep.names <- unique(dfilt$cells.info$culture)
#     #     logFC.df <- matrix(NA, nrow = nrow(dfilt), ncol = 2*length(rep.names) + 1)
#     #     colnames(logFC.df) <- c(paste0("logFC.", unique(dfilt$cells.info$culture)), "Ave.logFC")
#     #     treatment.names <- unique(dfilt$cells.info$perturbed)
#     #     for(i in 1:nrow(rp.results)) {
#     #         ## Code for one gene to do get logFC in each replicate experiment
#     #         dofit <- TRUE
#     # 
#     #         isexpr <- dfilt$isexpr_cpm[i,]
#     #         culture <- dfilt$cells.info$culture[isexpr]
#     #         treatment <- dfilt$cells.info$perturbed[isexpr]
#     #         if( (length(unique(culture)) > 1) & (length(unique(treatment)) > 1) ) {
#     #             design <- model.matrix(~ culture + treatment)
#     #             colnames(design) <- c("(Intercept)", "Expt", "Metformin")
#     #         }
#     #         else if( (length(unique(culture)) > 1) & (length(unique(treatment)) <= 1) ) {
#     #             warning("Gene expressed in only one treatment. Cannot do test on treatment.")
#     #             rp.results$N[i] <- sum(isexpr)
#     #             dofit <- FALSE
#     #         } else if( length(unique(treatment)) > 1 ) {
#     #             warning(paste0("Gene ", rownames(rp.results)[i], " expressed in only one culture/experiment (culture ", unique(culture), "). Results should be treated with caution."))
#     #             design <- model.matrix(~ treatment)
#     #             colnames(design) <- c("(Intercept)", "Metformin")
#     #         } else {
#     #             warning("Gene expressed in only one culture/experiment and one treatment. No DE results possible.")
#     #             rp.results$N[i] <- sum(isexpr)
#     #             dofit <- FALSE
#     #         }
#     #         if( dofit ) { ## Test to see if we should actually do the testing or not
#     #             x <- dfilt$log2cpm.pc0.25[i, isexpr]
#     #             rp.results$Ave.CPM[i] <- mean(x)
#     #             for( j in 1:length(rep.names) ) {
#     #                 ave.expr.1 <- mean(x[perturbed == treatment.names[1]] & culture == rep.names[j])
#     #                 ave.expr.2 <- mean(x[perturbed == treatment.names[2]] & culture == rep.names[j])
#     #                 logFC.df[, j] <- NA ### FIX ME
#     #             }
#     #         }
#     #     }
#     ## Compute average logFC and logFC-rankings for all genes
# }
# 
# #' Conduct RankProduct test for a data matrix
# #'
# #' @param data_matrix numeric matrix with expression values
# #' (log-scale) for features (rows) and samples/cells (columns).
# #' @param group character vector or factor defining the two groups to
# #' be compared
# #' @param block character vector or factor defining the experimental
# #' blocks within which the two groups are to be compared
# #' @param delta character string indication approximate p-value to be
# #' returned for rank-product results: "upper", "lower", or "geometric"
# #' (default; gives approximate p-value)
# #' @param compute_pvals logical, should p-values be computed? Default
# #' is \code{FALSE} as for experiments with many genes or cells,
# #' p-values will be extremely expensice or impossible to compute
# #' @return results object summarising RankProduct results
# #' @export
# testRankProduct.matrix <- function(data_matrix, group, block = NULL, delta = "geometric", compute_pvals = FALSE) {
#     ngenes <- nrow(data_matrix)
#     group <- as.factor(group)
#     block <- as.factor(block)
#     ## If no blocking, then assume one source for the data --> easy case
#     if(is.null(block) | nlevels(block) == 1) {
#         message("Treating study as having one experimental block.")
#         ## get logFCs
#         logFC <- getLogFCOneBlock(data_matrix, group)
#         ## Compute median logFC to report for each gene
#         median_logFC <- Biobase::rowMedians(logFC)
#         ## Get ranks for each pairwise comparison of cells across groups
#         ## Take log-ranks: true RP is likely to cause numeric overflow
#         lranks <- getRanksOneBlock(logFC, group, logged = TRUE)
#         lranks_up_grp1 <- lranks[[1]]
#         lranks_up_grp2 <- lranks[[2]]
#     }
#     else {
#         message("Computing rank-products across ", nlevels(block), " experimental blocks.")
#         blk_idx <- 1
#         for(blk in levels(block)) {
#             message("Block ", blk_idx, "...")
#             select_blk <- block == blk
#             this_grp <- group[select_blk]
#             ## get logFCs
#             this_logFC <- getLogFCOneBlock(data_matrix[, select_blk], this_grp)
#             ## Get ranks for each pairwise comparison of cells across groups
#             ## Take log-ranks: true RP is likely to cause numeric overflow
#             this_lranks <- getRanksOneBlock(this_logFC, this_grp, logged = TRUE)
#             if( blk_idx == 1 ) {
#                 logFC <- this_logFC
#                 lranks_up_grp1 <- this_lranks[[1]]
#                 lranks_up_grp2 <- this_lranks[[2]]
#             }
#             else {
#                 logFC <- cbind(logFC, this_logFC)
#                 lranks_up_grp1 <- cbind(lranks_up_grp1, this_lranks[[1]])
#                 lranks_up_grp2 <- cbind(lranks_up_grp2, this_lranks[[2]])
#             }
#             blk_idx <- blk_idx + 1
#         }
#     }
#     ## Compute median logFC to report for each gene
#     median_logFC <- Biobase::rowMedians(logFC)
#     ## Get mean log(rank-products) for up-regulation in group1 and up-regulation in group2
#     ## This is equivalent to the log of the geometric mean of the ranks
#     lrp_up_grp1 <- rowMeans(lranks_up_grp1)
#     lrp_up_grp2 <- rowMeans(lranks_up_grp2)
#     ## Report the smallest mean log-rank-product for each gene
#     lrp <- pmin(lrp_up_grp1, lrp_up_grp2)
#     ## Report whether expression higher in group1 or group2
#     label1 <- paste0(levels(group), collapse = ">")
#     label2 <- paste0(rev(levels(group)), collapse = ">")
#     direction <- rep(label1, length(lrp))
#     direction[lrp_up_grp1 > lrp_up_grp2] <- label2
#     ## Compute approximate p-values for geometric mean of rank-products, as if k = 1
#     if (compute_pvals) {
#         ## Compute the number of "replicates" if computing p-values
#         nreps <- ncol(logFC)
#         pvals <- calcRankProdPValsBounds(exp(lrp), ngenes, nreps, delta = delta)
#     }
#     else
#         pvals <- rep(NA, ngenes)
#     ## Define output
#     out_table <- data.frame(Median_logFC = median_logFC, Dir_of_Effect = direction,
#                             Geom_Mean_Rank = exp(lrp), P_Value = pvals)
#     rownames(out_table) <- rownames(data_matrix)
#     out_table
# }
# 
# 
# 
# 
