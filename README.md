
<!-- README.md is generated from README.Rmd. Please edit that file -->

# che

<!-- badges: start -->

<!-- badges: end -->

Convex Hull Ensembles for Species Distribution Modelling

A preprint exists at:

Hill, M., Caley, P., Camac, J., Elith, J., & Barry, S. (2022). A novel method accounting for predictor uncertainty and model transferability of invasive species distribution models. BioRxiv : The Preprint Server for Biology. https://doi.org/10.1101/2022.03.14.483865


## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mattecologist/che")
```

Note: this current version has mixed performance on different systems.
The `parallel` package is used to cycle through the pairwise models, and
I’ve only been testing it on Linux and Mac systems so far, so there may
be a need to test this on Windows. Tested on 4-core and 12-core
processors, the later gives much more stable performance.

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

![alt text here](man/figures/report_schematic.png) Schematic of new
ensemble modelling technique. A . Variable pairs of temperature and
precipitation using the bioclim variable are used in simple GAMs and
then assessed for relative model fit using AUC (and other metrics, see
text). Green squares are the pairs kept, red squares with crosses are
omitted. B . The variable pairs selected are used to create
two-dimensional environmental envelopes capturing all the presence
points (green points) inside all available environments (green + orange
points), using a minimum convex polygon (convex hulls (CHE)), or by use
of bounding boxes defined through the BIOCLIM algorithm (BBE). C . Each
of the two-variable envelopes from B are projected to a new geographical
surface (as a gridded raster) and then stacked (summed. The results are
then averaged, and the resulting surface is a continuous scale from 0-1
reflecting the level of certainty (closer to 1) that a given raster cell
falls in the environmental limits across the given predictor variables.

## Tutorial

`vignette("che")` is a short tutorial on how to use the CHE functions in
this package
