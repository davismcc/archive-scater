# Wrapper functions for the computation of p-values for exact rank-product tests
# Davis McCarthy, July 2014

library(matrixStats)

#' Conduct rank-product test
#'
#' @param dfilt DGEList object containing data and so on that has been
#' filtered and expressed cells determined
#' @param nmax.genes integer giving the maximum number of genes for
#' which to compute exact rank-product p-values
#' @param max.rprod.exact integer giving the maximum rank-product for
#' which to compute exact p-values. Exact p-values are very slow to
#' compute for large rank-products, and converge to approximate
#' p-values.
#' @return data frame with mean rank, rank product, number of ways of
#' obtaining given rank product, approximate p-value and exact p-value
#' @export
testRankProduct <- function(dfilt, nmax.genes = 500, max.rprod.exact = 100000) {
    navec <- rep(NA, nrow(dfilt))
    rp.results <- data.frame(Ave.CPM = navec, Ave.logFC = navec, N = navec)
    rownames(rp.results) <- rownames(dfilt)
    rep.names <- unique(dfilt$cells.info$culture)
    logFC.df <- matrix(NA, nrow = nrow(dfilt), ncol = 2*length(rep.names) + 1)
    colnames(logFC.df) <- c(paste0("logFC.", unique(dfilt$cells.info$culture)), "Ave.logFC")
    treatment.names <- unique(dfilt$cells.info$perturbed)
    for(i in 1:nrow(rp.results)) {
        ## Code for one gene to do get logFC in each replicate experiment
        dofit <- TRUE

        isexpr <- dfilt$isexpr_cpm[i,]
        culture <- dfilt$cells.info$culture[isexpr]
        treatment <- dfilt$cells.info$perturbed[isexpr]
        if( (length(unique(culture)) > 1) & (length(unique(treatment)) > 1) ) {
            design <- model.matrix(~ culture + treatment)
            colnames(design) <- c("(Intercept)", "Expt", "Metformin")
        }
        else if( (length(unique(culture)) > 1) & (length(unique(treatment)) <= 1) ) {
            warning("Gene expressed in only one treatment. Cannot do test on treatment.")
            rp.results$N[i] <- sum(isexpr)
            dofit <- FALSE
        } else if( length(unique(treatment)) > 1 ) {
            warning(paste0("Gene ", rownames(rp.results)[i], " expressed in only one culture/experiment (culture ", unique(culture), "). Results should be treated with caution."))
            design <- model.matrix(~ treatment)
            colnames(design) <- c("(Intercept)", "Metformin")
        } else {
            warning("Gene expressed in only one culture/experiment and one treatment. No DE results possible.")
            rp.results$N[i] <- sum(isexpr)
            dofit <- FALSE
        }
        if( dofit ) { ## Test to see if we should actually do the testing or not
            x <- dfilt$log2cpm.pc0.25[i, isexpr]
            rp.results$Ave.CPM[i] <- mean(x)
            for( j in 1:length(rep.names) ) {
                ave.expr.1 <- mean(x[perturbed == treatment.names[1]] & culture == rep.names[j])
                ave.expr.2 <- mean(x[perturbed == treatment.names[2]] & culture == rep.names[j])
                logFC.df[, j] <- NA ### FIX ME
            }
        }
    }
## Compute average logFC and logFC-rankings for all genes
}

#' Conduct RankProduct test for a data matrix
#'
#' @param data_matrix numeric matrix with expression values
#' (log-scale) for features (rows) and samples/cells (columns).
#' @param group character vector or factor defining the two groups to
#' be compared
#' @param block character vector or factor defining the experimental
#' blocks within which the two groups are to be compared
#' @param delta character string indication approximate p-value to be
#' returned for rank-product results: "upper", "lower", or "geometric"
#' (default; gives approximate p-value)
#' @return results object summarising RankProduct results
#' @export
testRankProduct.matrix <- function(data_matrix, group, block = NULL, delta = "geometric") {
    ngenes <- nrow(data_matrix)
    group <- as.factor(group)
    ## If no blocking, then assume one source for the data --> easy case
    if(is.null(block) | nlevels(block) == 1) {
        message("Treating study as having one experimental block.")
        ## get logFCs for 
        logFC <- getLogFCOneBlock(data_matrix, group)
        ## Compute median logFC to report for each gene
        median_logFC <- matrixStats::rowMedians(logFC)
        ## Get ranks for each pairwise comparison of cells across groups
        ## Take log-ranks: true RP is likely to cause numeric overflow
        lranks <- getRanksOneBlock(logFC, group, logged = TRUE)
        ## Get mean log(rank-products) for up-regulation in group1 and up-regulation in group2
        ## This is equivalent to the log of the geometric mean of the ranks
        lrp_up_grp1 <- rowMeans(lranks[[1]])
        lrp_up_grp2 <- rowMeans(lranks[[2]])
        ## Report the smallest mean log-rank-product for each gene
        lrp <- pmin(lrp_up_grp1, lrp_up_grp2)
        ## Report whether expression higher in group1 or group2
        direction <- rep(names(lranks)[1], length(lrp))
        direction[lrp_up_grp1 > lrp_up_grp2] <- names(lranks)[2]
        direction <- gsub("-", " > ", direction)
        ## Compute approximate p-values for geometric mean of rank-products, as if k = 1
        pvals <- calcRankProdBounds(exp(lrp), ngenes, 1, delta = delta)
        ## Define output
        out_table <- data.frame(Median_logFC = median_logFC, Dir_of_Effect = direction,
                                Geom_Mean_Rank = exp(lrp), P_Value = pvals)
    }
    else {
        message("Computing rank-products across ", nlevels(block), " experimental blocks.")
        stop("Multiple blocks not yet implemented.")
    }
    rownames(out_table) <- rownames(data_matrix)
    out_table
}


#' Get ranks from logFCs for one block
#'
#' Get ranks both for up-regulation in group1 and up-regulation in
#' group2
#'
#' @param logFC_matrix numeric matrix of logFCs computed as group1 -
#' group2
#' @param group character vector or factor defining the experimental
#' groups being compared
#' @param logged boolean: return raw ranks if \code{FALSE} or
#' log(ranks) if \code{TRUE}
#' @export
getRanksOneBlock <- function(logFC_matrix, group, logged = FALSE) {
    group <- as.factor(group)
    label1 <- paste0(levels(group), collapse = "-")
    label2 <- paste0(rev(levels(group)), collapse = "-")
    ## We assume that logFCs have been computed such that group1 = levels(group)[1] and logFC is group1 - group2
    ## Ranks in ascending order, so make logFC negative so that large positive logFCs end up at the top of the ranking
    rank_matrix_up_grp1 <- apply(-logFC_matrix, 2, rank)
    ## Rank by up-regulation in group2
    rank_matrix_up_grp2 <- apply(logFC_matrix, 2, rank)
    if(logged) {
        rank_matrix_up_grp1 <- log(rank_matrix_up_grp1)
        rank_matrix_up_grp2 <- log(rank_matrix_up_grp2)
    }
    ranks <- list(rank_matrix_up_grp1, rank_matrix_up_grp2)
    names(ranks) <- c(label1, label2)
    ranks
}


#' Get logFC values for replicates within one experimental block
#'
#' @param data_matrix numeric matrix with expression values
#' (log-scale) for features (rows) and samples/cells (columns).
#' @param group character vector or factor defining the two groups to
#' be compared
#' @return matrix of logFC values for all pairwise comparisons of
#' samples in group1 and group2
#' @export
getLogFCOneBlock <- function(data_matrix, group) {
    ## Check that group has the correct number of levels: only two groups currently supported
    group <- as.factor(group)
    if(nlevels(group) > 2)
        stop("More than two groups to be compared. Can only compare between two groups, currently.")
    ## Define indices for groups 1 and 2
    idx1 <- which(group == levels(group)[1])
    idx2 <- which(group == levels(group)[2])
    ngroup1 <- length(idx1)
    ngroup2 <- length(idx2)
    ## Define number of genes, comparisons and FC matrix
    ngenes <- nrow(data_matrix)
    ncells <- ncol(data_matrix)
    k <- ngroup1 * ngroup2
    logFC <- matrix(NA, nrow = ngenes, ncol = k)
    print("Hello!!")
    ## Get logFC for all pairwise comparisons of genes in group1 and group2
    for(i in seq_len(ngroup1)) {
        idx_this_insertion <- ((i - 1) * ngroup2 + 1):(i * ngroup2)
        insert_this <- data_matrix[, idx1[i]] - data_matrix[, idx2]
        logFC[, idx_this_insertion] <- insert_this
        ## colnames(logFC)[idx_this_insertion] <- paste(colnames(data_matrix)[idx1[i]],
        ##                                               colnames(data_matrix)[idx2], sep = "-")
        
    }
    ## Return logFC matrix
    logFC
}



#' Calculate exact rank-product p-values
#'
#' @param ranks a matrix of logFC-ranks, one row per gene, one column
#' per "experiment", i.e. pair of samples
#' @param nmax.genes integer giving the maximum number of genes for
#' which to compute exact rank-product p-values
#' @param max.rprod integer giving the maximum rank-product for
#' which to compute p-values. Exact p-values are very slow to
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
calcRankProdPVals <- function(ranks, nmax.genes = 500, max.rprod = 1e30) {
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
    ## for( i in 1:nmax.genes ) {
    ##     cat("Computing p-value for top gene", i, ",", names(rp.topgenes)[i], "\n")
    ##     this.rp <- rp.topgenes[i]
    ##     cat("Rank-product is:", this.rp, "\n")
    ##     if( this.rp <= max.rprod.exact ) {
    ##         while( n <= this.rp ) {
    ##             total <- total + piltzcount(n, nexp, ngenes)
    ##             n <- n + 1
    ##         }
    ##         counts[i] <- total
    ##     }
    ##     pvals.approx[i] <- calcRankProdPvalsApprox(this.rp, nexp, ngenes)
    ## }
    ## pvals <- counts / ngenes^nexp
    ## Unfortunately, more accurate approximate p-values take an impossibly long time to compute
    ## pvals.approx <- calcRankProdBounds(rp.topgenes, nexp, ngenes)
    pvals.approx <- calcRankProdPvalsApprox(rp.topgenes, nexp, ngenes)
    meanrank <- rp.topgenes^(1/nexp) # Get geometric mean of ranks
    out <- data.frame(Mean.Rank = meanrank, Rank.Product = rp.topgenes, Approx.P.Value = pvals.approx)
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


#' Calculate bounds on p-values for rank-product statistics
#'
#'
#' This function computes bounds on the p-value for rank products.
#'
#'
#' @param rho a vector of integers corresponding to the rank products
#' for which one wishes to compute the p-value.
#' @param n the number of features (e.g. genes).
#' @param k the number of replicates.
#' @param delta a character string indicating whether an upper bound
#' ('upper'), lower bound ('lower'), or the default,geometric
#' approximation (i.e. the geometric mean of the upper and lower
#' bounds) ('geometric') should be computed.
#'
#' @return a vector of p-values, one for each rank product.
#'
#' @details The exact p-value is guaranteed to be in between the lower
#' and the upper bound. The geometric mean of the two bounds can be
#' used as an approximation. Each bound is a piecewise continuous
#' function of the rank product. The different pieces each have an
#' analytic form, the parameters of which can be computed recursively.
#'
#' This implementation closely follows the description in Heskes,
#' Eisinga, Breitling: "A fast algorithm for determining bounds and
#' accurate approximate p-values of the rank product statistic for
#' replicate experiments", further referred to as HEB. More
#' specifically, this R function corresponds to the recursive variant,
#' sketched as pseudocode in the additional material of HEB.
#'
#' We thank HEB for making their R code available. Their original code
#' has been lightly modified for use in \pkg{scater}.
#'
#' @export
#' @examples
#' calcRankProdBounds(c(10, 100, 1000), 2000, 10, "upper")
calcRankProdBounds <- function(rho, n, k, delta = "geometric") {
    ## INPUT HANDLING
    if(any(rho > n^k) || any(rho < 1))
        stop('rho out of bounds')
    if(is.numeric(delta) == FALSE) {
        if(delta == 'geometric') {
            temp1 <- calcRankProdBounds(rho, n, k, 'upper')
            temp2 <- calcRankProdBounds(rho, n, k, 'lower')
            pvalue <- sqrt(temp1*temp2)   # geometric mean of upper and lower bound
            return(pvalue)
        }
        else {
            delta <- switch(delta,
                            upper = 1,        # for computing upper bound
                            lower = 0)        # for computing lower bound
        }
    }
    ## COMPUTE INTERVALS THAT CONTAIN THE RANK PRODUCTS
    logn <- log(n)
    allj <- ceiling(-(log(rho) / logn) + k)   # index specifying the interval that contains rho 
    minj <- min(allj)                     # lowest interval index
    maxj <- max(allj)                     # highest interval index
    ## INITIALIZE PARAMETERS
    param <- matrix(list(), nrow = k+1, ncol = maxj+1)
    for(i in 1:(k + 1)){
        for(j in 1:(maxj + 1)){
            param[[i, j]] <- list(a = c(), b = c(), c = c(), d = c(), e = c())
        }
    }
    ## param is a matrix of lists; each element of param is a list with values for the parameters
    ## a through e, which correspond to the parameters alpha through epsilon in HEB;
    ## specifially, param[[i+1,j+1]]$a corresponds to alpha_{i,j} in HEB, etc, where the offset
    ## of 1 is introduced to be able to represent, for example, alpha_{0,0};
    ## a, b, and c can be vectors (with possibly different lengths for different i and j),
    ## d and e are scalars
    ## COMPUTE PARAMETERS
    for(j in minj:maxj) {
        param <- updateParam(param, n, k, j, delta)
    }
    ## call to the function updateParam which recursively computes all parameters that are needed
    ## to calculate the p-value for a rank product rho that lies in the interval with index j
    ## COMPUTE RANK PRODUCTS GIVEN PARAMETERS
    k1 <- 1+k
    G <- rep(0, length(rho))   # G is a vector of the same length as rho,
                                        # for each rho bounding the number of rank products 
    for(j in minj:maxj) {
        j1 <- 1+j
        iii <- which(allj == j)         # indices of all rank products that fall in interval j:
                                    # bounds for these rank products can be computed with
                                    # the same set of parameters                                    
        thisrho <- rho[iii]
        thisparam <- param[[k1, j1]]
        thisG <- thisparam$e
        if(j != 0) {
            nrho <- length(thisrho)
            nterms <- length(thisparam$a)
            thisG <- thisG + thisparam$d * thisrho
            d1 <- matrix(thisparam$c) %*% thisrho
            d2 <- matrix(rep(log(thisrho), nterms), nrow = nterms, byrow = TRUE) -
                t(matrix(rep(logn * (k - j + thisparam$b), nrho), nrow = nrho, byrow = TRUE))
            d3 <- t(matrix(rep(thisparam$a, nrho), nrow = nrho, byrow = TRUE)) 
            thisG <- thisG + colSums(d1 * (d2^d3))
        }
                                        # the 10 lines above implement equation (8) in HEB
        G[iii] <- thisG
    }
    pvalue <- G/n^k
    return(pvalue)
}

#' Update parameters for computation of rank-product p-values
#'
#' This subroutine updates the current set of parameters to make sure
#' that the parameters corresponding to k replicates and the j'th
#' interval are included.
#'
#' @param param a matrix of lists, where each element of param is a
#' list with values for the parameters a through e; these parameters
#' specify the functional form of the bound; a, b, and c are all
#' vectors of unknown length, d and e are scalars.
#' @param n the number of features (e.g. genes).
#' @param k the number of replicates for which we need to compute the
#' corresponding parameters.
#' @param j the index of the interval for which we need to compute the
#' corresponding parameters.
#' @param delta 0 for the lower bound and 1 for the upper bound.
#'
#' @return A possibly updated set of parameters, at least including
#' those corresponding to (k,j).
#'
#' @details This subroutine makes sure that the parameters
#' corresponding to k replicates and a rank product within the j'th
#' interval are already present. If they already are (because
#' calculated before), it does not compute anything. Otherwise, it
#' recursively computes all parameters that are needed to arrive at
#' the parameters for (k,j).
#'
#' This implementation closely follows HEB, in particular equations (9) through (11).
#' 
updateParam <- function(param, n, k, j, delta) {
  k1 <- 1+k
  j1 <- 1+j
  if(length(param[[k1,j1]]$e) == 0) {  # apparently empty, so needs to be calculated
      if(j == 0) {   # initializing G_{k0}
          param[[k1,j1]]$e <- n^k
          param[[k1,j1]]$d <- 0
                                        # the 2 lines above implement equation (11) in HEB
      }
      else {
          k0 <- k1-1
          j0 <- j1-1
          param <- updateParam(param, n, k-1, j-1, delta)
                                        # checking that the parameters for (k-1,j-1) that are needed to compute the
                                        # parameters for (k,j) are indeed available; if not, they are themselves computed
          param00 = param[[k0, j0]]
          newa0 = param00$a + 1
          newb0 = param00$b
          newc0 = param00$c / newa0
          param11 = param00
                                        # the 5 lines above predefine some parameters common to equations (9) and (10) in HEB
          if(k == j){ # updates for G_{kk}
              param11$e <- (1 - delta) * (1 - param00$e)
              param11$d <- delta * param00$d + param00$e
              param11$a <- c(1, param00$a, newa0)
              param11$b <- c(0, param00$b, newb0)
              param11$c <- c(param00$d, delta * param00$c, newc0)
                                        # the 5 lines above implement equation (10) in HEB
          }
          else {  # updates for G_{kj}, j < k
              param <- updateParam(param, n, k-1, j, delta)
                                        # checking that the parameters for (k-1,j) that are needed to compute the
                                        # parameters for (k,j) are indeed available; if not, they are themselves computed
              param01 <- param[[k0, j1]]
              logn <- log(n)
              lognnkj <- (k - j) * logn
              newa1 <- param01$a + 1
              newa <- c(newa0, newa1)
              newb <- c(newb0, param01$b)
              newc <- c(newc0, -param01$c / newa1)
              param11$e <- n * param01$e + (delta - 1) * (param00$e - param01$e)
              lognminb <- c(-1 * param00$b * logn, (1 - param01$b) * logn)
              param11$d <- delta * param00$d + (1 - delta) * param01$d / n + 
                  (param00$e - param01$e) / exp(lognnkj) - sum(newc * (lognminb^newa))
              param11$a <- c(1, 1, param00$a, param01$a, newa)
              param11$b <- c(0, 1, param00$b, param01$b, newb)
              param11$c <- c(param00$d, -param01$d, delta * param00$c, (1-delta) * param01$c / n, newc)
                                        # the 15 lines above implement equation (9) in HEB
          }
          param[[k1, j1]] <- makeUnique(param11)
    # although not strictly necessary, the a, b and c vectors can possibly be shortened by
    # restricting oneselves to unique combinations of a and b values
      }
  }
  return(param)
}

#' Update the parameters so that they are unique
#'
#' This subroutine updates the parameters for a specific number of replicates and interval
#' such that it contains only unique combinations of the parameters a and b.
#'
#' @param param a single list with values for the parameters a through
#' e; these parameters specify the functional form of the bound; a, b,
#' and c are all vectors of unknown length, d and e are scalars.
#' 
#' @return A possibly updated and then more concise set of parameters
#' containing only unique combinations of the parameters a and b.
#'
#' @details While updating the vectors a and b, one may end up with
#' the exact same combinations of a and b. Given the functional form
#' of the bound, the representation can then be made more concise by
#' simply adding the corresponding elements of c.
makeUnique <- function(param) {
    ab <- t(rbind(param$a, param$b))
    uniqueab <- unique(ab)
    nunique <- dim(uniqueab)[1]
    param$a <- t(uniqueab[, 1])
    param$b <- t(uniqueab[, 2])
    newc <- rep(0, nunique)
    for(i in 1:nunique) {
        iii <- intersect(which(ab[, 1] == uniqueab[i, 1]), which(ab[, 2] == uniqueab[i, 2]))
        newc[i] <- sum(param$c[iii])  
    }
    param$c <- newc
    return(param)
}




