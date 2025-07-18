% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BetaBinomial.R
\name{smooth_bb}
\alias{smooth_bb}
\title{Perform Beta-Binomial Smoothing}
\usage{
smooth_bb(n_persons, n_items, freq, rmoment, nparm, rel)
}
\arguments{
\item{n_persons}{The total number of examinees in the sample.}

\item{n_items}{The number of items on the test.}

\item{freq}{A numeric vector of the observed raw score frequencies.}

\item{rmoment}{A numeric vector of the first four raw score moments.}

\item{nparm}{The number of parameters for the beta distribution (2 or 4).}

\item{rel}{The reliability of the test (e.g., KR-20), used to calculate
Lord's k for the compound binomial model. Set to 0 for the simple model.}
}
\value{
A list containing all results from the smoothing process, including
  the smoothed density, beta parameters, fitted moments, and chi-square statistics.
}
\description{
This is the main wrapper function that orchestrates the entire beta-binomial
smoothing process. It estimates true score moments, fits a beta distribution
to those moments, and then calculates the smoothed observed score distribution.
}
\author{
B. A. Hanson & R. L. Brennan (Original C code), Google's Gemini (R translation)
}
