% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap.R
\name{bootstrap_parametric_biv}
\alias{bootstrap_parametric_biv}
\title{Generate a Parametric Bootstrap Sample (Bivariate)}
\usage{
bootstrap_parametric_biv(n, smoothed_bfd)
}
\arguments{
\item{n}{The total sample size.}

\item{smoothed_bfd}{A matrix of smoothed bivariate probabilities.}
}
\value{
A matrix representing the new bootstrap bivariate frequency distribution.
}
\description{
Creates a parametric bootstrap sample from a smoothed bivariate distribution.
This is an R replacement for `Parametric_boot_biv`.
}
