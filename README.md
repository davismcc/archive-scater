# scater: single-cell analysis toolkit for expression with R

This package contains useful tools for the analysis of single-cell
gene expression data using the statistical software R. The package places an
emphasis on tools for quality control, visualisation and pre-processing of data
before further downstream analysis.

We hope that `scater` fills a useful niche between raw RNA-sequencing
count data and more focused downstream modelling tools such as 
[monocle](http://www.bioconductor.org/packages/release/bioc/html/monocle.html), 
[scLVM](http://github.com/PMBio/scLVM),
[SCDE](http://pklab.med.harvard.edu/scde/index.html), 
[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html), 
[limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) and 
so on.


## Installation
This package currently lives on GitHub, so I recommend using Hadley Wickham's 
`devtools` package to install `scater` directly from GitHub. If you don't have 
`devtools` installed, then install that from CRAN (as shown below) and then run
the call to install `scater`:

```{r }
install.packages("devtools")
devtools::install_github("davismcc/scater")
```

Using the most recent version of R is strongly recommended (R 3.2.0 at the time
of writing). Effort has been made to ensure the package works with R >3.0, but 
the package has not been tested with R <3.1.1.

There are several other packages from CRAN and Bioconductor that `scater` uses, 
so you will need to have these packages installed as well. The CRAN packages 
should install automatically when `scater` is installed, but you will need to 
install the Bioconductor packages manually.

Not all of the following are strictly necessary, but they enhance the 
functionality of `scater` and are good packages in their own right. The commands
below should help with package installations.

CRAN packages:

```{r}
install.packages(c("data.table", ggplot2", "ggthemes", "knitr", "magrittr", 
                 "matrixStats", "plyr", "reshape2", "rjson", "testthat"))
```

Bioconductor packages:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase", "BiocGenerics", "biomaRt", "edgeR", "limma", "monocle", 
        "rhdf5"))
```

You might like to install `dplyr` from Hadley Wickham's GitHub site (`dplyr` was
not available on CRAN for R 3.2.0 when this was written):

```{r }
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")
```

We plan to contribute `scater` to Bioconductor in the near future.


## Getting started

The best place to start is the [vignette](http://htmlpreview.github.io/?http://github.com/davismcc/scater/blob/master/vignettes/vignette.html).


## Acknowledgements and disclaimer

The package leans heavily on previously published work and packages, namely 
[edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) and 
[limma](http://bioconductor.org/packages/release/bioc/html/limma.html). The `SCESet` is heavily inspired by the `CellDataSet` class from [monocle](http://www.bioconductor.org/packages/release/bioc/html/monocle.html).


<!---
It also uses and extends code for an approximate rank-product test by [Heskes et al (2014)](http://dx.doi.org/10.1186/s12859-014-0367-1).
-->


The package is currently in an Alpha state, so is to be used with appropriate 
caution. Please do try it, though, and contact me with bug reports, feedback, 
questions and suggestions to improve the package.

Davis McCarthy, April 2015
