% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap.R
\name{bootstrap_bivariate}
\alias{bootstrap_bivariate}
\title{Generate a Non-Parametric Bootstrap Sample for a Bivariate Distribution}
\usage{
bootstrap_bivariate(n, bfd)
}
\arguments{
\item{n}{The total sample size.}

\item{bfd}{A matrix representing the original bivariate frequency distribution.}
}
\value{
A matrix representing the new bootstrap bivariate frequency distribution.
}
\description{
Creates a bootstrap sample from an observed bivariate frequency distribution.
This is an R replacement for the C function `Boot_BSTATS`.
}
