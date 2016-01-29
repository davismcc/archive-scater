## tests for plotting functions

context("test plotPCA, plotTSNE and plotDiffusionMap")

test_that("we can produce PCA plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    
    expect_output(class(plotPCA(example_sceset)), "ggplot")
})

test_that("we can produce t-SNE plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    
    expect_output(class(plotTSNE(example_sceset)), "ggplot")
})

test_that("we can produce Diffusion Map plots with different expression values",
          {
              data("sc_example_counts")
              data("sc_example_cell_info")
              pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
              example_sceset <- newSCESet(countData = sc_example_counts, 
                                          phenoData = pd)
              example_sceset <- calculateQCMetrics(example_sceset)
              
              expect_output(class(plotDiffusionMap(example_sceset)), "ggplot")
          })

context("test plotExpression")

test_that("we can produce expression plots with different expression values", {
    data("sc_example_counts")
    data("sc_example_cell_info")
    pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
    example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
    example_sceset <- calculateQCMetrics(example_sceset)
    
    expect_output(class(plotExpression(example_sceset, 1:4, "Cell_Cycle")), "ggplot")
    expect_output(class(plotExpression(example_sceset, 1:4, "Gene_0004")), "ggplot")
})