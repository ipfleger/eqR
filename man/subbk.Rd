% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CubicSpline.R
\name{subbk}
\alias{subbk}
\title{Solve an Upper Triangular System (Backward Substitution)}
\usage{
subbk(u, b)
}
\arguments{
\item{u}{An upper triangular numeric matrix.}

\item{b}{A numeric vector representing the right-hand side of the equation.}
}
\value{
A numeric vector `x` which is the solution to the system.
}
\description{
Solves the linear system `U %*% x = b` for `x`, where `U` is an upper
triangular matrix. This is an R equivalent of the `subbk` function.
}
\details{
This function is a wrapper for R's built-in `backsolve()` function,
which is optimized for this task. The C version modifies the `b` vector
in place, while this R function returns the solution vector `x`.
}
\examples{
U <- matrix(c(2, 6, -8, 0, 1, 5, 0, 0, 3), nrow = 3)
b <- c(4, 18, -2) # This is a different 'b' than the subfd example
x <- subbk(U, b)
print(x)
}
\author{
Jaehoon Seol (Original C code), Google's Gemini (R translation)
}
