% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ERutilities.R
\name{er_dfpmin}
\alias{er_dfpmin}
\title{Minimize a Function using the BFGS Algorithm}
\usage{
er_dfpmin(par, fn, gr = NULL, ...)
}
\arguments{
\item{par}{A numeric vector of initial values for the parameters.}

\item{fn}{The function to be minimized.}

\item{gr}{(Optional) The gradient of the function `fn`.}

\item{...}{Additional arguments to be passed to `fn` and `gr`.}
}
\value{
An object of class "optim" containing the optimization results.
}
\description{
Minimizes a function using the Broyden-Fletcher-Goldfarb-Shanno (BFGS)
algorithm. This is a wrapper for R's powerful `optim` function and replaces
the C library's `er_dfpmin` and its internal line search `er_lnsrch`.
}
