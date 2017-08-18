## Testing SCESet methods

test_that("newSCESet fails as expected", {
    expect_error(newSCESet(),
                 "Deprecated")

})
