% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ERutilities.R
\name{nscores}
\alias{nscores}
\title{Calculate the Number of Scores in a Sequence}
\usage{
nscores(max, min, inc)
}
\arguments{
\item{max}{A numeric value for the maximum score in the sequence.}

\item{min}{A numeric value for the minimum score in the sequence.}

\item{inc}{A numeric value for the increment between scores.}
}
\value{
An integer representing the total number of scores in the sequence.
}
\description{
This function calculates the total number of scores (or categories) in a
sequence defined by a minimum score, a maximum score, and a constant
increment. It is an R translation of the `nscores` function from the
"Equating Recipes" C code by R. L. Brennan.
}
\details{
This function relies on the `loc()` function to find the zero-based index
of the maximum score and then adds 1 to get the total count.

The original C code notes that for best results, the input values should be
specified with high precision (e.g., at least eight digits), especially
when the increment is a repeating decimal.
}
\examples{
# For a sequence of scores from 10 to 20 with an increment of 2,
# the scores are 10, 12, 14, 16, 18, 20. There are 6 scores in total.
number_of_scores <- nscores(max = 20, min = 10, inc = 2)
print(number_of_scores)
#> [1] 6
}
\seealso{
\code{\link{loc}}
}
\author{
R. L. Brennan (Original C code), Google's Gemini (R translation)
}
