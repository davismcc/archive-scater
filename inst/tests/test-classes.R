## Testing functions for classes used ##

context("tests on inputs")

test_that("tests for presence of cellData or counts",{
    expect_that(newSCESet(), 
                throws_error("Require at least one of cellData or countData object"))
})

# test_that("tests for dat variable",{
#     set.seed(12345)
#     grp <- as.factor(rep(c(0,1),each=15))
#     
#     dat <- matrix(0,nrow=100,ncol=30)
#     expect_that(deFunction(dat,grp),throws_error("some genes have zero variance; t-test won't work"))
# })

# context("test on outputs")
# 
# test_that("test p-values are numeric and non-zero",{
#     set.seed(12345)
#     grp <- as.factor(rep(c(0,1),each=15))
#     dat <- matrix(matrix(rnorm(100*30)),nrow=100,ncol=30)
#     
#     expect_that(deFunction(dat,grp),is_a("numeric"))
#     expect_that(all(deFunction(dat,grp) > 0),is_true())
# })


