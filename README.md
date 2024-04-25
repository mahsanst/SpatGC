# SpatGC <img src="man/figure/logo.png" align="right" width="25%"/>


[![R-CMD-check](https://github.com/mahsanst/SpatGC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mahsanst/SpatGC/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/SpatGC)](https://cran.r-project.org/package=SpatGC)
[![Downloads per
month](https://cranlogs.r-pkg.org/badges/SpatGC)](https://cran.r-project.org/package=SpatGC)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SpatGC)](https://cran.r-project.org/package=SpatGC)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)


The R package *SpatGC* Provides a collection of functions for preparing data and fitting Bayesian count spatial regression models, with a specific focus on the Gamma-Count (GC) model. The GC model is well-suited for modeling dispersed count data, including under-dispersed or over-dispersed counts, or counts with equivalent dispersion, using Integrated Nested Laplace Approximations (INLA). The package includes functions for generating data from the GC model, as well as spatially correlated versions of the model. See [Nadifar, Baghishani, Fallah (2023)](https://link.springer.com/article/10.1007/s13253-023-00550-5)


## Installation

You can install the **stable** version from

[CRAN](https://cran.r-project.org/package=SpatGC).

``` s
install.packages('SpatGC', dependencies = TRUE)
```

You can install the **development** version from
[Github](https://github.com/mahsanst/SpatGC)

``` s
# install.packages("remotes")
remotes::install_github("mahsanst/SpatGC")
``` 


## To cite package `SpatGC` in publications use:

Nadifar, M., Baghishani, H. (2024) *SpatGC:
Bayesian Modeling of Spatail Count Data*. R package version 0.1.0,
<https://cran.r-project.org/package=SpatGC>.

A BibTeX entry for LaTeX users is

@Manual{SpatGC, title = {SpatGC: Bayesian Modeling of Spatail Count Data},
author = {Mahsa Nadifar and Hossein Baghishani}, year = {2024}, note = {R package version 0.1.0}, url =
{<https://cran.r-project.org/package=SpatGC>} }

## License

This package is free and open-source software, licensed under GPL-3.

