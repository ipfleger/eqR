% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CubicSpline.R
\name{inverse_cubic_poly}
\alias{inverse_cubic_poly}
\title{Find the Inverse of a Cubic Polynomial}
\usage{
inverse_cubic_poly(yvalue, x_left, a, b, c, d)
}
\arguments{
\item{yvalue}{The y-value for which to find the corresponding x-value.}

\item{x_left}{The left boundary of the interval for the polynomial.}

\item{a, b, c, d}{The coefficients of the cubic polynomial.}
}
\value{
The x-value corresponding to `yvalue`.
}
\description{
Finds the x-value for a given y-value for a monotonic cubic polynomial,
effectively computing the inverse.
}
\details{
This function uses R's `uniroot` function to find the root of the equation
`f(x) - y = 0`, which is a robust and efficient method.
}
