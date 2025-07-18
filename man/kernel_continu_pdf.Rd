% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_Equate.R
\name{kernel_continu_pdf}
\alias{kernel_continu_pdf}
\title{Kernel Smoothed Continuous Probability Density Function (PDF)}
\usage{
kernel_continu_pdf(x, scores, rel_freq, hx)
}
\arguments{
\item{x}{A numeric value at which to evaluate the PDF.}

\item{scores}{A numeric vector of the discrete score points.}

\item{rel_freq}{A numeric vector of the relative frequencies for each score.}

\item{hx}{A numeric value for the bandwidth parameter.}
}
\value{
The value of the kernel smoothed PDF at `x`.
}
\description{
Computes the kernel smoothed continuous PDF for a given score `x`, based on
the methodology from von Davier, Holland, & Thayer (2004).
}
\author{
Tianyou Wang (Original C code), Google's Gemini (R translation)
}
