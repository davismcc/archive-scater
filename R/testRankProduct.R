# Wrapper functions for the computation of p-values for exact rank-product tests
# Davis McCarthy, July 2014

require(matrixStats)

#' Calculate exact rank-product p-values
#' 
#' @param ranks a matrix of logFC-ranks, one row per gene, one column per "experiment", i.e. pair of samples
#' @param nmax.genes integer giving the maximum number of genes for which to compute exact rank-product p-values
#' @return data frame with mean rank, rank product, number of ways of obtaining given rank product and exact p-value
calcRankProdPVals <- function(ranks, nmax.genes = 500) {
    ngenes <- nrow(ranks)
    nexp <- ncol(ranks)
    if( nmax.genes > ngenes )
        nmax.genes <- ngenes
    ## Calculate rank products
    rank.products <- matrixStats::rowProds(ranks)
    max.rp <- max(rank.products)
    cat("Maximum rank-product is: ", max.rp, "\n")
    print("Summary of rank-products:")
    print(summary(rank.products))
    ## Restrict computation to top nmax.genes genes
    cat("Computing p-values for top ", nmax.genes, "genes.\n")
    o <- order(rank.products)
    rp.topgenes <- rank.products[o][1:nmax.genes]
    cat("Max rank-product for top ", nmax.genes, "is ", max(rp.topgenes), ".\n")
    counts <- rep(NA, nmax.genes)
    prev.rp <- 0
    total <- 0
    for( i in 1:nmax.genes ) {
        this.rp <- rp.topgenes[i]
        for( n in (prev.rp + 1):this.rp ) {
            total <- total + piltzcount(n, nexp, ngenes)
        }
        counts[i] <- prev.counts <- total  
    }
    pvals <- counts / ngenes^nexp
    meanrank <- rp.topgenes^(1/nexp) # Get geometric mean of ranks
    out <- data.frame(Mean.Rank = meanrank, Rank.Product = rp.topgenes, N.Ways = counts, P.Value = pvals)
    rownames(out) <- names(rp.topgenes)
    out
}
