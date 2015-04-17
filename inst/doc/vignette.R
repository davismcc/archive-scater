## ----knitr-options, echo=FALSE-------------------------------------------
library(knitr)
opts_chunk$set(fig.align='center', fig.width=6, fig.height=5)

## ----quickstart-load-data------------------------------------------------
library(scater)
data("sc_example_counts")
data("sc_example_cell_info")
ls()

## ----quickstart-make-sceset----------------------------------------------
pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
example_sceset

## ----filter-no-exprs-----------------------------------------------------
keep_gene <- rowSums(isExprs(example_sceset)) > 0
example_sceset <- example_sceset[keep_gene,]

## ----sceset-make-sceset-counts-only--------------------------------------
pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
gene_df <- data.frame(Gene = rownames(sc_example_counts))
rownames(gene_df) <- gene_df$Gene
fd <- new("AnnotatedDataFrame", data = gene_df)
example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd,
                            featureData = fd)
example_sceset

## ----sceset-make-sceset-exprs-only---------------------------------------
example2 <- newSCESet(cellData = edgeR::cpm.default(sc_example_counts))
pData(example2)
fData(example2)

## ----counts-accessor-----------------------------------------------------
counts(example2)[1:3, 1:6]

## ----exprs-accessor------------------------------------------------------
exprs(example2)[1:3, 1:6]

## ----isexprs-accessor----------------------------------------------------
isExprs(example2)[1:3, 1:6]

## ----sceset-define-is-exprs----------------------------------------------
isExprs(example2) <- calcIsExprs(example2, lowerDetectionLimit = 0.5, 
                               expr_data = "exprs")
isExprs(example2)[1:3, 1:6]

## ----sceset-add-count-data-----------------------------------------------
counts(example2) <- sc_example_counts
example2
counts(example2)[1:3, 1:6]

## ----sceset-demo-replacement---------------------------------------------
gene_df <- data.frame(Gene = rownames(sc_example_counts))
rownames(gene_df) <- gene_df$Gene
fd <- new("AnnotatedDataFrame", data = gene_df)
## replace featureData
fData(example_sceset) <- fd
## replace phenotype data
pData(example_sceset) <- pd
## replace expression data to be used
exprs(example_sceset) <- edgeR::cpm.default(counts(example_sceset), 
                                            prior.count = 5, log = TRUE)

## ----plot-expression-----------------------------------------------------
plotExpression(rownames(example_sceset)[1:6], example_sceset, 
               aes(x = Mutation_Status, y = log2_evals, colour = IsExpressed))

## ----plot-expression-theme-bw--------------------------------------------
theme_set(theme_bw())
plotExpression(rownames(example_sceset)[7:12], example_sceset, 
               aes(x = Mutation_Status, y = log2_evals, colour = IsExpressed),
                   show_median = TRUE, show_violin = TRUE, 
               xlab="Mutation Status")


## ----calc-qc-metrics-----------------------------------------------------
names(pData(example_sceset))
example_sceset <- calculateQCMetrics(example_sceset, control_genes=NULL)
names(pData(example_sceset))

## ----plot-qc-------------------------------------------------------------
keep_gene <- rowSums(isExprs(example_sceset)) > 0
example_sceset <- example_sceset[keep_gene,]
## Plot QC
plotQC(example_sceset, type="most-expressed")
plotQC(example_sceset, type="find-pcs", variable="coverage")
plotQC(example_sceset, type="expl", variable=c("coverage", "depth"))

## ----plot-pdata----------------------------------------------------------
plotPhenoData(example_sceset, aes(x=depth, y=coverage, colour=Mutation_Status))

## ----plot-fdata----------------------------------------------------------
plotFeatureData(example_sceset, aes(x=n_cells_exprs, y=prop_total_reads))

## ----plot-pdata-move-legend----------------------------------------------
plotPhenoData(example_sceset, aes(x=depth, y=coverage, colour=Mutation_Status)) +
    theme(legend.position="top")

## ----plot-pdata-move-legend-add-trend------------------------------------
plotPhenoData(example_sceset, aes(x=depth, y=coverage, colour=Mutation_Status)) +
    theme(legend.position="top") +
    stat_smooth(method="lm", se=FALSE, size=2, fullrange=TRUE)

