# scater: single-cell analysis tools for expression with R

This package contains useful tools for the analysis of single-cell
gene expression data using the statistical software R. The package places an
emphasis on tools for quality control, visualisation and pre-processing of data
before further downstream analysis.

## Installation
This package currently lives on GitHub, so I recommend using Hadley Wickham's 
`devtools` package to install `scater` directly from GitHub. If you don't have 
`devtools` installed, then install that from CRAN (as shown below) and then run
the call to install `scater`:

```{r }
install.packages("devtools")
devtools::install_github("davismcc/scater")
```

## Getting started

The best place to start is the [vignette](http://github.com/davismcc/scater/blob/master/vignettes/vignette.Rmd).

## Acknowledgements and disclaimer

The package leans heavily on previously published work and packages, namely [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) and [limma](http://bioconductor.org/packages/release/bioc/html/limma.html). It also uses and extends code for an approximate rank-product test by [Heskes et al (2014)](http://dx.doi.org/10.1186/s12859-014-0367-1).

The package is currently in an Alpha state, so is to be used with caution. 
Please do try it, though, and contact me with bug reports, feedback, questions
and suggestions to improve the package.

Davis McCarthy, April 2015
