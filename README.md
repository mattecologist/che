
<!-- README.md is generated from README.Rmd. Please edit that file -->

# che

<!-- badges: start -->

<!-- badges: end -->

Convex Hull Ensembles for Species Distribution Modelling

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mattecologist/che")
```

## Background

Species distribution models (SDMs) have become increasingly popular to
predict the potential geographic range of invasive species, by
correlating the known distribution of a species with environmental
variables, and transferring these relationships into novel environments
and geographic space. The SDM literature has largely focused on
improving the statistical fit of models when projecting into novel
geographic space. Specifically, by understanding how background size and
extent may influence model transferability, and how best to select
absence, pseudo-absence or background points for modelling. These
investigations are coupled with considerable examination of how model
evaluation methods and test scores are influenced by these choices. An
area that remains challenging and under-discussed is methods to choose
predictive variables that will result in models that transfer between
regions. This is challenging as non-causal variables may still show
significant associations with distributions.

Here we introduce methods to finesse this problem by using multiple
simple models to search for the most appropriate variables from a given
set, and then apply them as two variable envelopes in an ensemble
approach. This is in the `che` function of this package. We also include
modified approaches using `bioclim` variables and `range bagging`.

## Example

There will soon be a vignette demonstrating use of this package
