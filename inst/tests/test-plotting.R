## tests for plotting functions

context("test plotPCA and plotTSNE")

test_that("we can produce PCA plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    
    expect_output(class(plotPCA(example_sceset)), "ggplot")
})