# Compare and aggregate results from different ranking methods using exact rank-product tests
# Davis McCarthy, July 2014

library(matrixStats)

#' Conduct rank-product test
#'
#' @param results a matrix or data.frame object where columns provide
#' results values (e.g. p-values) to be used for ranking features
#' (e.g. genes) and ordering of results.
#' @param nmax.features integer giving the maximum number of features
#' (e.g. genes) for which to compute exact rank-product p-values
#' @param max.rprod.exact integer giving the maximum rank-product for
#' which to compute exact p-values. Exact p-values are very slow to
#' compute for large rank-products, and converge to approximate
#' p-values.
#' @return data frame with mean rank, rank product, number of ways of
#' obtaining given rank product, approximate p-value and exact p-value
#' @export
aggregateResults <- function(results, nmax_genes = 500, max_rankprod_exact = 100000) {
    ranks <- matrix(NA, nrow = nrow(results), ncol = ncol(results), dimnames = list(rownames(results), colnames(results)))
    ## Define ranking for each feature, for each method
    for(i in seq_len(ncol(results))) {
        ranks[, i] <- rank(results[, i], ties.method = "average")
    }
    rank_products <- matrixStats::rowProds(ranks)
    names(rank_products) <- rownames(ranks)
    ranks <- ranks[order(rank_products), ]
    aggregated_results <- calcRankProdPVals(ranks, nmax_genes, max_rankprod_exact)
    aggregated_results <- cbind(ranks[1:nmax_genes, ], aggregated_results)
    aggregated_results
}
    
#' Calculate exact rank-product p-values
#'
#' @param ranks a matrix of logFC-ranks, one row per gene, one column
#' per "experiment", i.e. pair of samples
#' @param nmax.genes integer giving the maximum number of genes for
#' which to compute exact rank-product p-values
#' @param max.rprod.exact integer giving the maximum rank-product for
#' which to compute exact p-values. Exact p-values are very slow to
#' compute for large rank-products, and converge to approximate
#' p-values.
#' @return data frame with mean rank, rank product, number of ways of
#' obtaining given rank product, approximate p-value and exact p-value
#' @export
#' @examples
#' ptm <- proc.time()
#' ranks <- replicate(3, sample(50))
#' rownames(ranks) <- paste0("gene.", 1:50)
#' colnames(ranks) <- paste0("exp.", 1:3)
#' results <- calcRankProdPVals(ranks)
#' results
#' proc.time() - ptm
calcRankProdPVals <- function(ranks, nmax.genes = 500, max.rprod.exact = 100000) {
    ngenes <- nrow(ranks)
    nexp <- ncol(ranks)
    if( nmax.genes > ngenes )
        nmax.genes <- ngenes
    ## Calculate rank products
    rank.products <- matrixStats::rowProds(ranks)
    names(rank.products) <- rownames(ranks)
    max.rp <- max(rank.products)
    cat("Maximum rank-product is: ", max.rp, "\n")
    print("Summary of rank-products:")
    print(summary(rank.products))
    ## Restrict computation to top nmax.genes genes
    cat("Computing p-values for top ", nmax.genes, "genes.\n")
    o <- order(rank.products)
    rp.topgenes <- rank.products[o][1:nmax.genes]
    cat("Max rank-product for top", nmax.genes, "genes is", max(rp.topgenes), ".\n")
    counts <- pvals.approx <- rep(NA, nmax.genes)
    prev.rp <- 0
    total <- 0
    n <- 1
    for( i in 1:nmax.genes ) {
        cat("Computing p-value for top gene", i, ",", names(rp.topgenes)[i], "\n")
        this.rp <- rp.topgenes[i]
        cat("Rank-product is:", this.rp, "\n")
        if( this.rp <= max.rprod.exact ) {
            while( n <= this.rp ) {
                total <- total + piltzcount(n, nexp, ngenes)
                n <- n + 1
            }
            counts[i] <- total
        }
        pvals.approx[i] <- calcRankProdPvalsApprox(this.rp, nexp, ngenes)
    }
    pvals <- counts / ngenes^nexp
    meanrank <- rp.topgenes^(1/nexp) # Get geometric mean of ranks
    out <- data.frame(Mean.Rank = meanrank, Rank.Product = rp.topgenes, N.Ways = counts, Approx.P.Value = pvals.approx, Exact.P.Value = pvals)
    rownames(out) <- names(rp.topgenes)
    out
}

#' Right-tailed gamma approximation to rank-product p-value
#'
#' @param r integer rank-product
#' @param k integer number of experiments
#' @param n integer number of genes
#' @return float aproximate p-value
#' @export
#' @examples
#' calcRankProdPvalsApprox(9720, 5, 500)
calcRankProdPvalsApprox <- function(r, k, n) {
    1 - pgamma(-log(r/(n+1)^k), k, scale = 1)
}

