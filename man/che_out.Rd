% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/che.R
\name{che_out}
\alias{che_out}
\title{Reconcile CHE output}
\usage{
che_out(model.out = model.out, AUC_threshold = 0.25)
}
\arguments{
\item{model.out}{The model output}

\item{AUC_threshold}{The AUC threshold to use to select layers to include}
}
\value{
A list of a single raster, the AUC scores, and the AUC value corresponding to the chosen threshold
}
\description{
Format and summarise the CHE model outputs
}
