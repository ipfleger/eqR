% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print_cll_cg}
\alias{print_cll_cg}
\title{Print a Summary of Continuized Log-Linear Equating for CG Design}
\usage{
print_cll_cg(title, eq, pdata, conversion_table)
}
\arguments{
\item{title}{A character string for the main title of the output.}

\item{pdata}{A list object (like the one returned by `smooth_bll`) containing
the setup and parameter information for the analysis. It should include
details about the smoothing model (e.g., polynomial degrees).}

\item{conversion_table}{A list or data frame containing the equated scores. It is
expected to have a column for the original raw scores and additional columns
for the equated scores from each method (e.g., "FE", "ChainedE").}
}
\value{
This function does not return a value. It prints its output to the
  console.
}
\description{
This function takes the results from a continuized log-linear equating run
for a Common-Item Non-Equivalent Groups (CG) design and prints a formatted
summary to the console. It is an R equivalent of the C library's `Print_CC`
function.
}
\examples{
# This is a conceptual example of how you would use it.
# Assume 'pdata_list' and 'results_df' are returned from your main wrapper.

# print_cll_cg(
#   title = "My Equating Run: Form B to Form A",
#   pdata = pdata_list,
#   conversion_table = results_df
# )
}
