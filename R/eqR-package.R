## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom stats approx approxfun cov dnorm median optim optimize quantile
#'   rmultinom sd setNames spline uniroot update
#' @importFrom utils head tail
#' @importFrom graphics abline legend lines plot polygon
#' @importFrom grDevices adjustcolor hcl.colors
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib eqR, .registration = TRUE
## usethis namespace: end
NULL

# Non-standard-evaluation symbols (formula terms in add_plan, a lazy default arg
# in a CLL statistic helper) that R CMD check cannot see bound.
utils::globalVariables(c("from", "to", "score_params"))
