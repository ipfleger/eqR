#' Perform Cholesky Decomposition
#'
#' @description
#' This function computes the Cholesky decomposition of a symmetric,
#' positive-definite matrix. It serves as an R equivalent to the `dpdch`
#' function from the "Equating Recipes" C library.
#'
#' @details
#' The function uses R's built-in `chol()` function, which is highly optimized
#' for this purpose. The original C function `dpdch` overwrites the input
#' matrix, whereas this R function returns a new matrix containing the upper
#' triangular factor of the decomposition.
#'
#' @param matx A symmetric, positive-definite numeric matrix.
#'
#' @return An upper triangular matrix `R` such that `t(R) %*% R` is equal to
#'   the original matrix `matx`.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' \dontrun{
#' A <- matrix(
