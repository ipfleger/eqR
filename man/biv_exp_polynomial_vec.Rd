% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CLL_Equate.R
\name{biv_exp_polynomial_vec}
\alias{biv_exp_polynomial_vec}
\title{Vectorized Bivariate Exponential Polynomial Calculation (Corrected)}
\usage{
biv_exp_polynomial_vec(bivar, x, y)
}
\arguments{
\item{bivar}{A `bivar` object from the smoothing process.}

\item{x}{A vector of scores for the first variable (u).}

\item{y}{A vector of scores for the second variable (v).}
}
\value{
A numeric vector of results from the exponential polynomial calculation.
}
\description{
A corrected, vectorized version that uses pre-calculated indices from the
bivar object to correctly and efficiently build the design matrix.
}
