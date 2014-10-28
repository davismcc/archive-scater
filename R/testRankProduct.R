## # Wrapper functions for the computation of p-values for exact rank-product tests
## # Davis McCarthy, July 2014

## library(matrixStats)

## #' Conduct rank-product test
## #'
## #' @param dfilt DGEList object containing data and so on that has been
## #' filtered and expressed cells determined
## #' @param nmax.genes integer giving the maximum number of genes for
## #' which to compute exact rank-product p-values
## #' @param max.rprod.exact integer giving the maximum rank-product for
## #' which to compute exact p-values. Exact p-values are very slow to
## #' compute for large rank-products, and converge to approximate
## #' p-values.
## #' @return data frame with mean rank, rank product, number of ways of
## #' obtaining given rank product, approximate p-value and exact p-value
## #' @export
## testRankProduct <- function(dfilt, nmax.genes = 500, max.rprod.exact = 100000) {
##     navec <- rep(NA, nrow(dfilt))
##     rp.results <- data.frame(Ave.CPM = navec, Ave.logFC = navec, N = navec)
##     rownames(rp.results) <- rownames(dfilt)
##     rep.names <- unique(dfilt$cells.info$culture)
##     logFC.df <- matrix(NA, nrow = nrow(dfilt), ncol = 2*length(rep.names) + 1)
##     colnames(logFC.df) <- c(paste0("logFC.", unique(dfilt$cells.info$culture)), "Ave.logFC")
##     treatment.names <- unique(dfilt$cells.info$perturbed)
##     for(i in 1:nrow(rp.results)) {
##         ## Code for one gene to do get logFC in each replicate experiment
##         dofit <- TRUE

##         isexpr <- dfilt$isexpr_cpm[i,]
##         culture <- dfilt$cells.info$culture[isexpr]
##         treatment <- dfilt$cells.info$perturbed[isexpr]
##         if( (length(unique(culture)) > 1) & (length(unique(treatment)) > 1) ) {
##             design <- model.matrix(~ culture + treatment)
##             colnames(design) <- c("(Intercept)", "Expt", "Metformin")
##         }
##         else if( (length(unique(culture)) > 1) & (length(unique(treatment)) <= 1) ) {
##             warning("Gene expressed in only one treatment. Cannot do test on treatment.")
##             rp.results$N[i] <- sum(isexpr)
##             dofit <- FALSE
##         } else if( length(unique(treatment)) > 1 ) {
##             warning(paste0("Gene ", rownames(rp.results)[i], " expressed in only one culture/experiment (culture ", unique(culture), "). Results should be treated with caution."))
##             design <- model.matrix(~ treatment)
##             colnames(design) <- c("(Intercept)", "Metformin")
##         } else {
##             warning("Gene expressed in only one culture/experiment and one treatment. No DE results possible.")
##             rp.results$N[i] <- sum(isexpr)
##             dofit <- FALSE
##         }
##         if( dofit ) { ## Test to see if we should actually do the testing or not
##             x <- dfilt$log2cpm.pc0.25[i, isexpr]
##             rp.results$Ave.CPM[i] <- mean(x)
##             for( j in 1:length(rep.names) ) {
##                 ave.expr.1 <- mean(x[perturbed == treatment.names[1]] & culture == rep.names[j])
##                 ave.expr.2 <- mean(x[perturbed == treatment.names[2]] & culture == rep.names[j])
##                 logFC.df[, j] <- NA ### FIX ME
##             }
##         }
##     }
## ## Compute average logFC and logFC-rankings for all genes
## }


## #' Get logFC values for replicates within on experiment
## getLogFC.oneexperiment <- function() {

## }



## #' Calculate exact rank-product p-values
## #'
## #' @param ranks a matrix of logFC-ranks, one row per gene, one column
## #' per "experiment", i.e. pair of samples
## #' @param nmax.genes integer giving the maximum number of genes for
## #' which to compute exact rank-product p-values
## #' @param max.rprod.exact integer giving the maximum rank-product for
## #' which to compute exact p-values. Exact p-values are very slow to
## #' compute for large rank-products, and converge to approximate
## #' p-values.
## #' @return data frame with mean rank, rank product, number of ways of
## #' obtaining given rank product, approximate p-value and exact p-value
## #' @export
## #' @examples
## #' ptm <- proc.time()
## #' ranks <- replicate(3, sample(50))
## #' rownames(ranks) <- paste0("gene.", 1:50)
## #' colnames(ranks) <- paste0("exp.", 1:3)
## #' results <- calcRankProdPVals(ranks)
## #' results
## #' proc.time() - ptm
## calcRankProdPVals <- function(ranks, nmax.genes = 500, max.rprod.exact = 100000) {
##     ngenes <- nrow(ranks)
##     nexp <- ncol(ranks)
##     if( nmax.genes > ngenes )
##         nmax.genes <- ngenes
##     ## Calculate rank products
##     rank.products <- matrixStats::rowProds(ranks)
##     names(rank.products) <- rownames(ranks)
##     max.rp <- max(rank.products)
##     cat("Maximum rank-product is: ", max.rp, "\n")
##     print("Summary of rank-products:")
##     print(summary(rank.products))
##     ## Restrict computation to top nmax.genes genes
##     cat("Computing p-values for top ", nmax.genes, "genes.\n")
##     o <- order(rank.products)
##     rp.topgenes <- rank.products[o][1:nmax.genes]
##     cat("Max rank-product for top", nmax.genes, "genes is", max(rp.topgenes), ".\n")
##     counts <- pvals.approx <- rep(NA, nmax.genes)
##     prev.rp <- 0
##     total <- 0
##     n <- 1
##     for( i in 1:nmax.genes ) {
##         cat("Computing p-value for top gene", i, ",", names(rp.topgenes)[i], "\n")
##         this.rp <- rp.topgenes[i]
##         cat("Rank-product is:", this.rp, "\n")
##         if( this.rp <= max.rprod.exact ) {
##             while( n <= this.rp ) {
##                 total <- total + piltzcount(n, nexp, ngenes)
##                 n <- n + 1
##             }
##             counts[i] <- total
##         }
##         pvals.approx[i] <- calcRankProdPvalsApprox(this.rp, nexp, ngenes)
##     }
##     pvals <- counts / ngenes^nexp
##     meanrank <- rp.topgenes^(1/nexp) # Get geometric mean of ranks
##     out <- data.frame(Mean.Rank = meanrank, Rank.Product = rp.topgenes, N.Ways = counts, Approx.P.Value = pvals.approx, Exact.P.Value = pvals)
##     rownames(out) <- names(rp.topgenes)
##     out
## }

## #' Right-tailed gamma approximation to rank-product p-value
## #'
## #' @param r integer rank-product
## #' @param k integer number of experiments
## #' @param n integer number of genes
## #' @return float aproximate p-value
## #' @export
## #' @examples
## #' calcRankProdPvalsApprox(9720, 5, 500)
## calcRankProdPvalsApprox <- function(r, k, n) {
##     1 - pgamma(-log(r/(n+1)^k), k, scale = 1)
## }

