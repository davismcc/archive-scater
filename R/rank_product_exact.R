## #   The exact probability distribution of the rank product statistics for
## #   replicated experiments, by Eisinga, Breitling & Heskes, FEBS Letters,
## #   January 2013

## # Modified by Davis McCarthy
## # April 2014


## #PILTZCOUNT Number of ways to construct a rank product
## #   piltzcount(r,k,n) returns the number of ways rank product r can be
## #   constructed from k experiments if n is the number of genes (and thus
## #   the maximum rank in a single experiment)
## #

## #---------- load CRAN packages -------------------------------------------

## #require(gmp)
## require(mgcv)

## #-------------------------------------------------------------------------

## #' Get the unique rows of a matrix
## #' 
## #' Save the unique rows of input matrix A in a new matrix C and create
## #' an index vector ia, containing the sequence numbers of the rows in A
## #' that were saved in C.
## #' 
## #' @param A matrix for which we obtain unique rows
## #' @return A list containing the matrix C, containing only the unique rows of A, and an index vector \code{ia} giving the indices of the unique rows of A in C.
## #' @keywords manip
## #' @export
## #' @examples
## #' A <- matrix(rep(1:4, each = 4), byrow = TRUE, ncol = 2)
## #' uniquerows(A)
## uniquerows <- function(A) {
##     C <- mgcv::uniquecombs(A) # uniquecombs() from package mgcv
##     ia <- vector("integer", length=nrow(C))
##     for( i in 1:nrow(C) ) {
##         for( j in 1:nrow(A) ) {
##             #score <- ifelse(C[i,] == A[j,], yes = 1, no = 0)
##             #if( sum(score) == ncol(A) )
##             if( identical(C[i,], A[j,]) ) {
##                 ia[i] <- j
##                 break
##             }
##         }
##     }
##     return(list(C, ia))
## } # end function uniquerows

## #------------------------------------------------------------------------

## #' Extract prime factors and powers of rank product.
## #' 
## #' Given an integer (in this context a rank product), factorize it into its prime factors and their powers.
## #' 
## #' @param n integer to be factorized
## #' @return data frame in which each row gives a unique prime factor and its power
## #' @keywords arithmetic
## #' @export
## #' @examples
## #' tidyfactor(281137500)
## ## tidyfactor <- function(n) {
## ##     if( n==1 )
## ##         return(c(1,0))
## ##     ## Use factorize() from package gmp.
## ##     pf <- as.integer(gmp::factorize(n))
## ##     f  <- unique(pf)
## ##     ## Compute multiplicities.
## ##     nprimes <- length(f)
## ##     m <- matrix(nrow=nprimes, ncol=1, 0)
## ##     for( i in 1:nprimes )
## ##         m[i] <- sum(pf == f[i])
## ##     return(cbind(f, m))
## ## } # end tidyfactor.

## #------------------------------------------------------------------------

## #' Obtain all divisors of a set of prime numbers and their powers
## #' 
## #' Obatin all divisors of $f_1^m_1 \times f_2^m_2 \times f_3^m_3 \ldots$, where $f_1, f_2, f_3 \ldots$ are primes and $m_1, m_2, m_3 \ldots$ are prime powers.
## #' 
## #' @param f.and.m matrix in which the first column gives the primes and the second column gives the prime powers, e.g. the output of \code{tidyfactor}
## #' @return sorted vector of divisors of $f_1^m_1 \times f_2^m_2 \times f_3^m_3 \ldots$
## #' @keywords arithmetic
## #' @export
## #' @seealso tidyfactor
## #' @examples
## #' f.and.m <- tidyfactor(281137500)
## #' tidydivisor(f.and.m)
## tidydivisor <- function(f.and.m) {
##     f <- f.and.m[, 1]
##     m <- f.and.m[, 2]
##     nprimes <- length(f)
##     d <- matrix(f[1]^(0:m[1]), ncol=1)
##     if( nprimes > 1 ) {
##           for( i in 2:nprimes ) {
##               d <- d %*% f[i]^(0:m[i])
##               d <- matrix(d, ncol=1)
##           }
##       }
##     d <- sort(d)
##     return(d)
## } # end tidydivisor

## #-------------------------------------------------------------------------

## #' Calculate binomial coefficient
## #' 
## #' Calculate the binomial coefficient.
## #' 
## #' @param n integer giving total number of items
## #' @param k integer giving number of items to choose
## #' @return integer value for the binomial coefficient
## #' @export
## #' @seealso choose
## #' @examples
## #' bincoeff(10, 2)
## bincoeff <- function(n, k) {
##     y <- exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1))
##     y <- round(y)
##     return(y)
## } # end bincoeff

## #-------------------------------------------------------------------------

## #' Calculate multinomial coefficient
## #' 
## #' Calculate multinomial coefficient
## #' 
## #' @param n integer giving total number of items
## #' @param x vector of number of items in each of the multinomial classes 
## #' @return integer value for the multinomial coefficient
## #' @export
## #' @seealso choose
## #' @examples
## #' multcoeff(10, c(2, 1, 3))
## multcoeff <- function(n, x) {
##     y <- exp(lgamma(n+1) - sum(lgamma(x+1)))
##     y <- round(y)
##     return(y)
## } # end multcoeff

## #-------------------------------------------------------------------------

## #' Number of ways to construct a rank-product
## #' 
## #' Compute the number of k-tuples of rank-product value r, given n genes, which is equal to the number of ways a rank-product r can be constructed from k experiments if n is the number of genes (and thus the maximum rank in a single experiment)
## #' 
## #' @param r integer rank-product value
## #' @param k integer number of experiments
## #' @param n integer number of genes
## #' @return integer value for the number of ways to construct the rank-product
## #' @export
## #' @examples
## #' #--- Example 1.
## #' # Calculate for n=500 the number of (k=)5-tuples with rankproduct rp=9720.
## #' piltzcount(9720, 5, 500)
## #' 
## #' # Gamma distribution approximation of p value for rank product rp=9720.
## #' righttailgamma <- function(r, k, n) 1 - pgamma(-log(r/(n+1)^k), k, scale=1)
## #' 
## #' righttailgamma(9720, 5, 500)
## #' #--- Example 2.
## #' # Calculate for n=500 the number of (k=)5-tuples with rankproduct rp<=9720.
## #' # Display the exact p-value.
## #' total <- 0
## #' for( n in 1:9720 ) {
## #'     if( n %% 1000 == 0 )
## #'         cat("n=", n, "\n")
## #'     total <- total + piltzcount(n, 5, 500)
## #' }
## #' print(total)
## #' total / 500^5
## piltzcount <- function(r, k, n) {

##     #cat ("r=", r, " k=", k, " n=", n, "\n")

##     ## Check for trivial cases.
##     if( k < 1 ) {
##         h <- 0
##         return(h)
##     } 
##     if( r == 1 | r > n^k ) {
##         h <- ifelse(r == 1, 1, 0)
##         return(h)
##     }
##     if( k == 1 ) {
##         h <- 1
##         return(h)
##     }

##     ## Compute prime factorization.
##     fm <- tidyfactor(r)
##     f <- fm[,1]; m <- fm[,2]
##     nprimes <- length(f)

##     # Calculate first term: number of ordered k-tuples with product r.
##     h <- 1
##     for( t in 1:nprimes )
##         h <- h * bincoeff(m[t]+k-1, k-1)

##     #cat("First term of h: ", h, "\n")

##     if( r > n ) {

##         ## Get all divisors of r.
##         d <- tidydivisor(fm)
##         ## Select divisors larger than n.
##         bigones <- d[d>n]  # bigones is a row vector.
##         nbig <- length(bigones)
##         #print("bigones")
##         #print(bigones)

##         ## Construct possible combinations of divisors.
##         divset <- list()
##         bbb <- list()
##         divset[[1]] <- 1
##         bbb[[1]] <- matrix(0, nrow=1, ncol=nbig)

##         ## Compute maximum number of divisors that can be divided out: smax.
##         smax <- ceiling(log(r)/log(n)) - 1

##         for( s in 1:smax ) {

##             #print("s"); print(s)
##             newdivset <- matrix(divset[[s]], ncol=1) %*% bigones
##             #print("newdivset")
##             #print(newdivset)

##             ij <- as.matrix(which(r - newdivset*floor(r/newdivset) < 0.5, arr.ind=TRUE))
##             #print(ij)
##             i <- ij[,1]
##             j <- ij[,2]
##             nnew <- length(i)
##             #cat("nnew=", nnew,"\n")

##             if( nnew > 0 ) {

##                 divset[[s+1]] <- matrix(0, ncol=nnew)
##                 bbb[[s+1]]    <- matrix(0, nrow=nnew, ncol=nbig)

##                 for( t in 1:nnew ) {
##                     #print("***")
##                     #print(divset[[s]][i[t]])
##                     divset[[s+1]][t] <- divset[[s]][i[t]] * bigones[j[t]]
##                     bbb[[s+1]][t,] <- bbb[[s]][i[t],]
##                     bbb[[s+1]][t,j[t]] <- bbb[[s+1]][t,j[t]] + 1
##                 }
##                 #print("divset[[s+1]]")
##                 #print(divset[[s+1]])
##                 #print("bbb")
##                 #print(bbb)

##                 temp <- uniquerows(bbb[[s+1]])
##                 bbb[[s+1]] <- temp[[1]];  iii <- temp[[2]]
##                 #print("iii")
##                 #print(iii)
##                 divset[[s+1]] <- divset[[s+1]][iii]
##                 #print("divset")
##                 #print(divset)

##                 for( t in 1:length(iii) ) {
##                     hextra <- piltzcount(r / divset[[s+1]][t], k-s, r)
##                     #cat("hextra=", hextra, "s= ", s, "divset[[s+1]][t]= ", divset[[s+1]][t], "\n")
##                     add_or_subt <- bincoeff(k,s) * multcoeff(s,bbb[[s+1]][t,]) * hextra
##                     #cat("add_or_subt= ", add_or_subt, "\n")
##                     ## Update h, number of ways of getting rank-product
##                     h <- h + (-1)^s * add_or_subt
##                 }

##             } else {
##                 break
##             }

##         } # end s in 1:smax

##     } #end if (r > n)

##     return(h)

## } # end piltzcount

## #-------------------------------------------------------------------------





