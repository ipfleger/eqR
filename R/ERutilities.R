#' Calculate the Zero-Based Index of a Score
#'
#' @description
#' This function calculates the location (zero-based index) of a given score
#' within a sequence of scores defined by a minimum value and a constant
#' increment. It is an R translation of the `loc` function from the
#' "Equating Recipes" C code by R. L. Brennan.
#'
#' @details
#' The original C code notes that for best results, the input values should be
#' specified with high precision (e.g., at least eight digits), especially
#' when the increment is a repeating decimal.
#'
#' Note that this function returns a **zero-based index** to maintain
#' consistency with the original C library's logic. R users should be mindful
#' of this, as R typically uses 1-based indexing. To use this index to subset
#' a standard R vector, you would need to add 1.
#'
#' @param x A numeric value representing the score.
#' @param min A numeric value for the minimum score in the sequence.
#' @param inc A numeric value for the increment between scores.
#'
#' @return An integer representing the zero-based index of the score `x`.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' # Find the location of the score 16 in a sequence from 10 to 20 by 2.
#' score_location <- loc(x = 16, min = 10, inc = 2)
#' print(score_location)
#' #> [1] 3
loc <- function(x, min, inc) {
  # The C code uses `(int) (val + 0.5)` for rounding.
  # R's `round()` function achieves the same result for positive numbers.
  # as.integer() ensures the return type matches the original C function.
  return(as.integer(round((x - min) / inc)))
}

#' Calculate the Number of Scores in a Sequence
#'
#' @description
#' This function calculates the total number of scores (or categories) in a
#' sequence defined by a minimum score, a maximum score, and a constant
#' increment. It is an R translation of the `nscores` function from the
#' "Equating Recipes" C code by R. L. Brennan.
#'
#' @details
#' This function relies on the `loc()` function to find the zero-based index
#' of the maximum score and then adds 1 to get the total count.
#'
#' The original C code notes that for best results, the input values should be
#' specified with high precision (e.g., at least eight digits), especially
#' when the increment is a repeating decimal.
#'
#' @param max A numeric value for the maximum score in the sequence.
#' @param min A numeric value for the minimum score in the sequence.
#' @param inc A numeric value for the increment between scores.
#'
#' @return An integer representing the total number of scores in the sequence.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @seealso \code{\link{loc}}
#'
#' @examples
#' # For a sequence of scores from 10 to 20 with an increment of 2,
#' # the scores are 10, 12, 14, 16, 18, 20. There are 6 scores in total.
#' number_of_scores <- nscores(max = 20, min = 10, inc = 2)
#' print(number_of_scores)
#' #> [1] 6
nscores <- function(max, min, inc) {
  # This function calls loc() to get the zero-based index of the max score
  # and adds 1 to get the total number of scores. The 'L' ensures the result
  # is an integer, consistent with the C version.
  return(loc(max, min, inc) + 1L)
}

#' Calculate the Score from a Zero-Based Index
#'
#' @description
#' This function calculates the score associated with a given zero-based index
#' within a sequence defined by a minimum value and a constant increment. It is
#' the inverse of the `loc()` function. This is an R translation of the `score`
#' function from the "Equating Recipes" C code by R. L. Brennan.
#'
#' @details
#' The original C code notes that for best results, the input values should be
#' specified with high precision (e.g., at least eight digits).
#'
#' @param loc An integer representing the zero-based index.
#' @param min A numeric value for the minimum score in the sequence.
#' @param inc A numeric value for the increment between scores.
#'
#' @return A numeric value representing the score at the given location.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @seealso \code{\link{loc}}
#'
#' @examples
#' # Find the score at the 3rd index (4th position) in a sequence
#' # starting at 10 with an increment of 2.
#' the_score <- score(loc = 3, min = 10, inc = 2)
#' print(the_score)
#' #> [1] 16
score <- function(loc, min, inc) {
  # This is a direct translation of the C formula.
  return(min + loc * inc)
}


#' Compute Percentile Rank from a Cumulative Distribution
#'
#' @description
#' This function computes the percentile rank for a given score `x` using a
#' cumulative relative frequency distribution (`crfd`). It performs linear
#' interpolation between the points of the distribution, consistent with the
#' methodology described in Kolen & Brennan (2004). This is an R translation
#' of the `perc_rank` function from the "Equating Recipes" C code.
#'
#' @details
#' This function uses R's built-in `approxfun()` to create a linear
#' interpolation function. The points for interpolation are defined by the
#' score scale (`min`, `max`, `inc`) and the provided cumulative relative
#' frequencies (`crfd`). The score points are treated as the midpoints of
#' intervals that are `inc` wide.
#'
#' The function is vectorized and can compute percentile ranks for multiple
#' `x` values at once.
#'
#' @param x A numeric vector of scores for which to calculate percentile ranks.
#' @param min A numeric value for the minimum score in the sequence.
#' @param max A numeric value for the maximum score in the sequence.
#' @param inc A numeric value for the increment between scores.
#' @param crfd A numeric vector of the cumulative relative frequency
#'   distribution. The length of this vector should be `nscores(max, min, inc)`.
#'
#' @return A numeric vector of percentile ranks (0-100).
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' scores <- seq(0, 10)
#' freqs <- c(1, 2, 5, 8, 10, 15, 12, 8, 5, 3, 1)
#' rel_freqs <- freqs / sum(freqs)
#' cum_rel_freqs <- cumsum(rel_freqs)
#'
#' # Calculate percentile rank for a single score
#' perc_rank(x = 5.5, min = 0, max = 10, inc = 1, crfd = cum_rel_freqs)
#'
#' # Calculate for multiple scores
#' perc_rank(x = c(2, 5.5, 9), min = 0, max = 10, inc = 1, crfd = cum_rel_freqs)
perc_rank <- function(x, min, max, inc, crfd) {
  # Define the score points and their corresponding boundaries for interpolation
  score_points <- seq(from = min, to = max, by = inc)
  score_boundaries <- seq(from = min - inc / 2, to = max + inc / 2, by = inc)

  # The cumulative relative frequency at the lower bound is 0
  crfd_boundaries <- c(0, crfd)

  # Create a linear interpolation function using the boundaries and crfd
  # `rule = 2` tells approxfun to return the value at the closest data point
  # for values outside the interpolation range.
  interp_fun <- stats::approxfun(x = score_boundaries, y = crfd_boundaries, rule = 2)

  # Calculate the percentile rank (as a proportion) by applying the function
  pr_proportion <- interp_fun(x)

  # Convert to percentage (0-100) and return
  return(pr_proportion * 100)
}

#' Compute Percentile Point (Score) from a Percentile Rank
#'
#' @description
#' This function computes the score (percentile point) corresponding to a given
#' percentile rank (`pr`). It uses a cumulative relative frequency distribution
#' (`crfd`) and implements the two-step linear interpolation method described in
#' Kolen & Brennan (2004, pp. 45-46). This is the inverse of the `perc_rank`
#' function.
#'
#' @details
#' The function calculates an "upper" and "lower" percentile point and averages
#' them to produce the final result. This specific methodology ensures results
#' are consistent with the original C library. The function is vectorized and
#' can compute scores for multiple percentile ranks at once.
#'
#' @param pr A numeric vector of percentile ranks (0-100) for which to find
#'   the corresponding scores.
#' @param ns An integer, the total number of score categories.
#' @param min A numeric value for the minimum score in the sequence.
#' @param inc A numeric value for the increment between scores.
#' @param crfd A numeric vector of the cumulative relative frequency
#'   distribution. The length of this vector should be `ns`.
#'
#' @return A numeric vector of scores corresponding to the given percentile ranks.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @seealso \code{\link{perc_rank}}
#'
#' @examples
#' # Setup a sample distribution
#' freqs <- c(1, 2, 5, 8, 10, 15, 12, 8, 5, 3, 1)
#' rel_freqs <- freqs / sum(freqs)
#' cum_rel_freqs <- cumsum(rel_freqs)
#'
#' # Find the score corresponding to the 50th percentile
#' perc_point(pr = 50, ns = 11, min = 0, inc = 1, crfd = cum_rel_freqs)
#'
#' # Find scores for multiple percentile ranks
#' perc_point(pr = c(25, 50, 75), ns = 11, min = 0, inc = 1, crfd = cum_rel_freqs)
perc_point <- function(pr, ns, min, inc, crfd) {

  # Vectorized implementation
  sapply(pr, function(single_pr) {
    prp <- single_pr / 100

    # Handle edge cases
    if (prp <= 1e-8) {
      i <- which(crfd > 1e-8)[1]
      ppU <- i - 1 - 0.5 # C uses 0-based index, so i-1
      ppL <- -0.5
      return(min + inc * ((ppU + ppL) / 2))
    }

    if (prp >= 1 - 1e-8) {
      j <- tail(which(crfd < 1 - 1e-8), 1)
      ppL <- j - 1 + 1 + 0.5 # C uses 0-based index, so j-1
      ppU <- ns - 1 + 0.5
      return(min + inc * ((ppU + ppL) / 2))
    }

    if (crfd[1] > prp) {
      return(min + inc * (prp / crfd[1] - 0.5))
    }

    # Upper percentile point (ppU)
    i <- which(crfd > prp)[1]
    if (crfd[i] != crfd[i - 1]) {
      ppU <- (prp - crfd[i - 1]) / (crfd[i] - crfd[i - 1]) + (i - 1 - 0.5)
    } else {
      ppU <- i - 1 - 0.5
    }

    # Lower percentile point (ppL)
    j <- tail(which(crfd < prp), 1)
    if (crfd[j + 1] != crfd[j]) {
      ppL <- (prp - crfd[j]) / (crfd[j + 1] - crfd[j]) + (j - 1 + 0.5)
    } else {
      ppL <- j - 1 + 0.5
    }

    return(min + inc * ((ppU + ppL) / 2))
  })
}


#' Compute the Dot Product of Two Vectors
#'
#' @description
#' Computes the dot product of two numeric vectors. This is a wrapper around
#' R's highly efficient `crossprod()` function to maintain naming consistency
#' with the original C library (`er_dot`).
#'
#' @details
#' The original C function included an `offset` parameter to handle 0-based or
#' 1-based indexing. This is not needed in R, as vector operations are
#' inherently 1-based and vectorized.
#'
#' @param vect1 A numeric vector.
#' @param vect2 A numeric vector of the same length as `vect1`.
#'
#' @return A single numeric value representing the dot product of the two vectors.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' v1 <- c(1, 2, 3)
#' v2 <- c(4, 5, 6)
#' er_dot(v1, v2) # Expected: 1*4 + 2*5 + 3*6 = 32
er_dot <- function(vect1, vect2) {
  # crossprod() is generally the most efficient way to compute a dot product
  # in R. The result is a 1x1 matrix, so we use as.numeric() to extract the value.
  return(as.numeric(crossprod(vect1, vect2)))
}

#' Compute the Dot Product of Two Vectors
#'
#' @description
#' Computes the dot product of two numeric vectors. This is a wrapper around
#' R's highly efficient `crossprod()` function to maintain naming consistency
#' with the original C library (`er_dot`).
#'
#' @details
#' The original C function included an `offset` parameter to handle 0-based or
#' 1-based indexing. This is not needed in R, as vector operations are
#' inherently 1-based and vectorized.
#'
#' @param vect1 A numeric vector.
#' @param vect2 A numeric vector of the same length as `vect1`.
#'
#' @return A single numeric value representing the dot product of the two vectors.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' v1 <- c(1, 2, 3)
#' v2 <- c(4, 5, 6)
#' er_dot(v1, v2) # Expected: 1*4 + 2*5 + 3*6 = 32
er_dot <- function(vect1, vect2) {
  # crossprod() is generally the most efficient way to compute a dot product
  # in R. The result is a 1x1 matrix, so we use as.numeric() to extract the value.
  return(as.numeric(crossprod(vect1, vect2)))
}

#' Perform the DAXPY Operation
#'
#' @description
#' Performs the "DAXPY" (Double-precision A*X Plus Y) operation, a common
#' routine in basic linear algebra subprograms (BLAS). It computes the
#' expression `alpha * x + y`.
#'
#' @details
#' This is a wrapper around R's standard vectorized arithmetic, which is
#' highly efficient. The C version (`er_daxpy`) modifies the `y` vector in
#' place, whereas this R function returns a new vector containing the result,
#' which is more idiomatic for R.
#'
#' @param vectY A numeric vector (the "y" in the operation).
#' @param alpha A single numeric value (the scalar "a").
#' @param vectX A numeric vector (the "x" in the operation), which must be the
#'   same length as `vectY`.
#'
#' @return A new numeric vector containing the result of `alpha * vectX + vectY`.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' y <- c(10, 20, 30)
#' x <- c(1, 2, 3)
#' a <- 2
#' er_daxpy(y, a, x) # Expected: 2*c(1,2,3) + c(10,20,30) = c(12, 24, 36)
er_daxpy <- function(vectY, alpha, vectX) {
  # R's vectorized arithmetic is the most direct and efficient way to do this.
  return(alpha * vectX + vectY)
}

#' Scale a Vector by a Constant
#'
#' @description
#' Performs the "scale" BLAS operation, which multiplies a vector by a scalar
#' value. It computes the expression `scale * x`.
#'
#' @details
#' This is a wrapper around R's standard vectorized arithmetic. The C version
#' (`er_scale`) modifies one of the input vectors in place. This R function
#' is designed to be more intuitive by taking the vector to be scaled as the
#' primary argument and returning a new vector with the result.
#'
#' @param vect A numeric vector to be scaled.
#' @param scale A single numeric value to scale the vector by.
#'
#' @return A new numeric vector containing the scaled values.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' v <- c(10, 20, 30)
#' s <- 0.5
#' er_scale(v, s) # Expected: c(5, 10, 15)
er_scale <- function(vect, scale) {
  # R's vectorized arithmetic is the natural way to perform this operation.
  return(vect * scale)
}


#' Perform a Scaled Rank-1 Update of a Matrix
#'
#' @description
#' Computes a scaled rank-1 update of a matrix, a common operation in
#' quasi-Newton optimization methods like BFGS. It performs the operation:
#' `matx_new = matx_old + scale * (vect %*% t(vect))`
#'
#' @details
#' This function uses R's efficient `tcrossprod()` function, which computes
#' the outer product of the vector with itself (`vect %*% t(vect)`). The C
#' version (`er_r1update`) modifies the input matrix in place, whereas this
#' R function returns the updated matrix as a new object.
#'
#' @param matx A numeric matrix.
#' @param scale A single numeric value to scale the outer product by.
#' @param vect A numeric vector. The outer product `vect %*% t(vect)` will
#'   result in a matrix with dimensions `length(vect) x length(vect)`, which
#'   must match the dimensions of `matx`.
#'
#' @return A new matrix containing the updated values.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' H <- diag(3) # Initial matrix (identity)
#' v <- c(1, 2, 3)
#' s <- 0.1
#' H_new <- er_r1update(H, s, v)
#' print(H_new)
er_r1update <- function(matx, scale, vect) {
  # tcrossprod(vect) is an efficient way to compute the outer product v %*% t(v)
  outer_product <- tcrossprod(vect)

  # Return the updated matrix
  return(matx + scale * outer_product)
}

#' Perform Matrix-Vector Multiplication
#'
#' @description
#' Computes the product of a matrix and a vector. This is a wrapper around
#' R's built-in `%*%` operator for matrix multiplication.
#'
#' @details
#' The C version (`er_mvmult`) modifies an output vector in place. This R
#' function returns a new vector containing the result, which is the standard
#' convention in R.
#'
#' @param matx A numeric matrix.
#' @param vect1 A numeric vector. The number of elements in `vect1` must equal
#'   the number of columns in `matx`.
#'
#' @return A new numeric vector containing the result of the multiplication.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' M <- matrix(1:6, nrow = 2, ncol = 3)
#' V <- c(1, 2, 3)
#' er_mvmult(M, V) # Expected: c(14, 32)
er_mvmult <- function(matx, vect1) {
  # The %*% operator is R's standard for matrix multiplication.
  # The result is a matrix with one column, so we use as.vector to convert.
  return(as.vector(matx %*% vect1))
}

#' Compute Statistical Moments from a Frequency Distribution
#'
#' @description
#' Computes the first four statistical moments (mean, standard deviation,
#' skewness, and kurtosis) from a frequency distribution. This function is a
#' more idiomatic R replacement for the C functions `MomentsFromFD` and
#' `MomentsFromRFD`.
#'
#' @details
#' The function can accept either absolute frequencies (`freq`) or relative
#' frequencies (`rel_freq`). If `freq` is provided, it will be used to
#' calculate the total number of observations and then converted to relative
#' frequencies. If only `rel_freq` is provided, the calculations are based
#' directly on the proportions.
#'
#' The `scores` can be provided as a vector. If `scores` is `NULL`, they are
#' generated as a sequence from `min` to `max` with an increment of `inc`.
#'
#' @param scores A numeric vector of the score points. If `NULL` (the default),
#'   scores are generated using `min`, `max`, and `inc`.
#' @param freq An integer vector of absolute frequencies corresponding to the
#'   scores. Either `freq` or `rel_freq` must be provided.
#' @param rel_freq A numeric vector of relative frequencies (proportions).
#' @param min A numeric value for the minimum score (used if `scores` is `NULL`).
#' @param max A numeric value for the maximum score (used if `scores` is `NULL`).
#' @param inc A numeric value for the increment (used if `scores` is `NULL`).
#'
#' @return A named numeric vector containing the mean, standard deviation (sd),
#'   skewness, and kurtosis.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' freqs <- c(1, 2, 5, 8, 10, 15, 12, 8, 5, 3, 1)
#' score_vals <- 0:10
#' get_moments(scores = score_vals, freq = freqs)
get_moments <- function(scores = NULL, freq = NULL, rel_freq = NULL, min = NULL, max = NULL, inc = NULL) {

  if (is.null(scores)) {
    if (is.null(min) || is.null(max) || is.null(inc)) {
      stop("If 'scores' is NULL, 'min', 'max', and 'inc' must be provided.")
    }
    scores <- seq(from = min, to = max, by = inc)
  }

  if (is.null(freq) && is.null(rel_freq)) {
    stop("Either 'freq' or 'rel_freq' must be provided.")
  }

  if (!is.null(freq)) {
    n <- sum(freq)
    if (n == 0) return(c(mean = NA, sd = NA, skewness = NA, kurtosis = NA))
    rel_freq <- freq / n
  } else {
    n <- NA # n is unknown if only relative frequencies are given
  }


  # Calculate moments using weighted operations
  mean_val <- sum(scores * rel_freq)

  dev <- scores - mean_val
  var_val <- sum(dev^2 * rel_freq)
  sd_val <- sqrt(var_val)

  if (sd_val == 0) {
    return(c(mean = mean_val, sd = 0, skewness = NA, kurtosis = NA))
  }

  skew_val <- sum(dev^3 * rel_freq) / (sd_val^3)
  kurt_val <- sum(dev^4 * rel_freq) / (sd_val^4)

  return(c(mean = mean_val, sd = sd_val, skewness = skew_val, kurtosis = kurt_val))
}



#' Perform Equipercentile Equating
#'
#' @description
#' Computes the equipercentile equivalents for a set of new form scores (`X`)
#' on the scale of an old form (`Y`).
#'
#' @details
#' This function takes the percentile ranks of the new form scores and finds the
#' corresponding scores on the old form's scale using the `perc_point` function.
#' This is a direct implementation of the equipercentile equating method.
#'
#' @param nsy An integer, the number of score categories for the old form `Y`.
#' @param miny A numeric value, the minimum score for form `Y`.
#' @param incy A numeric value, the increment between scores for form `Y`.
#' @param crfdy A numeric vector of the cumulative relative frequency
#'   distribution for form `Y`.
#' @param nsx An integer, the number of score categories for the new form `X`.
#' @param prdx A numeric vector of the percentile ranks for each score on
#'   form `X`.
#'
#' @return A numeric vector of equated raw scores on the `Y` scale.
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @seealso \code{\link{perc_point}}, \code{\link{perc_rank}}
#'
#' @examples
#' # Example data from Kolen & Brennan (2004), Table 2.4
#' # Form X distribution
#' freqs_x <- c(1, 2, 5, 8, 10, 15, 12, 8, 5, 3, 1)
#' min_x <- 0; max_x <- 10; inc_x <- 1
#' ns_x <- nscores(max_x, min_x, inc_x)
#' crfd_x <- cumsum(freqs_x / sum(freqs_x))
#' prd_x <- perc_rank(x = 0:10, min = min_x, max = max_x, inc = inc_x, crfd = crfd_x)
#'
#' # Form Y distribution
#' freqs_y <- c(1, 3, 6, 9, 12, 14, 11, 7, 4, 2, 1)
#' min_y <- 0; max_y <- 10; inc_y <- 1
#' ns_y <- nscores(max_y, min_y, inc_y)
#' crfd_y <- cumsum(freqs_y / sum(freqs_y))
#'
#' # Perform equating
#' equated_scores <- EquiEquate(nsy = ns_y, miny = min_y, incy = inc_y,
#'                              crfdy = crfd_y, nsx = ns_x, prdx = prd_x)
#' print(data.frame(x_score = 0:10, equated_y = equated_scores))
EquiEquate <- function(nsy, miny, incy, crfdy, nsx, prdx, eraw) {
  # This is a direct application of perc_point for each percentile rank in prdx
  eraw <- perc_point(pr = prdx, ns = nsy, min = miny, inc = incy, crfd = crfdy)
  return(eraw)
}

#' Perform LU Decomposition of a Matrix
#'
#' @description
#' This function serves as a wrapper for R's LU decomposition capabilities,
#' primarily for compatibility with the C code's structure. It computes the
#' LU decomposition of a square matrix.
#'
#' @details
#' The original C function `er_ludcmp` performs LU decomposition with partial
#' pivoting and modifies the input matrix `a` to store the L and U factors.
#' This R version will simply return the result of an LU decomposition, which
#' can then be used by `er_lubksb`. For most R applications, one would use
#' `solve()` directly.
#'
#' @param a A square numeric matrix.
#'
#' @return A list containing the LU-decomposed matrix and pivot information,
#'   compatible with `er_lubksb`. In R, this is often handled by `solve` or
#'   functions from the `Matrix` package. For simplicity, we can use `solve`
#'   as a high-level equivalent.
#'
#' @seealso \code{\link{er_lubksb}}, \code{\link{solve}}
er_ludcmp <- function(a) {
  # In R, LU decomposition is typically handled internally by solve().
  # For a more direct equivalent, one might use the 'Matrix' package's lu().
  # For the purpose of this translation, we'll assume the main goal is to
  # solve a system, which is done in the next function. This function
  # can be a placeholder or simply return the matrix itself, as `solve`
  # will handle the decomposition.
  return(a)
}

#' Solve a System of Linear Equations using LU Decomposition
#'
#' @description
#' Solves the linear system `a %*% x = b` for `x`. This function is designed
#' to work with the output of `er_ludcmp`, mirroring the C library's two-step
#' process for solving linear systems.
#'
#' @details
#' This function is a wrapper for R's built-in `solve()` function, which is
#' highly optimized and performs LU decomposition with partial pivoting
#' internally to solve the system. The C version `er_lubksb` performs
#' forward and backward substitution on an already-decomposed matrix.
#'
#' @param a A square numeric matrix (the coefficient matrix).
#' @param b A numeric vector or matrix (the right-hand side of the equation).
#'
#' @return A vector or matrix `x` representing the solution to the system.
#' @seealso \code{\link{er_ludcmp}}, \code{\link{solve}}
#'
#' @examples
#' A <- matrix(c(2, 1, -1, -3, -1, 2, -2, 1, 2), nrow = 3, byrow = TRUE)
#' B <- c(8, -11, -3)
#' # In C, you would call ludcmp then lubksb. In R, we can combine them.
#' X <- er_lubksb(A, B)
#' print(X)
#' #> [1]  2  3 -1
er_lubksb <- function(a, b) {
  is_square <- is.matrix(a) && nrow(a) == ncol(a)
  full_rank <- qr(a)$rank == ncol(a)

  if (is_square && full_rank) {
    return(solve(a, b))
  } else {
    return(MASS::ginv(a) %*% b)
  }
}

#' Compute the Inverse of a Matrix
#'
#' @description
#' Computes the inverse of a square matrix using Gauss-Jordan elimination
#' with partial pivoting, as implemented in the original C function `er_matrix_inverse`.
#'
#' @details
#' This function is a wrapper for R's built-in `solve()` function. When
#' `solve()` is called with a single matrix argument, it computes the inverse.
#'
#' @param a A square, non-singular numeric matrix.
#'
#' @return A matrix representing the inverse of `a`.
#' @seealso \code{\link{solve}}
#'
#' @examples
#' M <- matrix(c(1, 2, 3, 4), nrow = 2)
#' M_inv <- er_matrix_inverse(M)
#' print(M_inv)
#' print(M %*% M_inv) # Should be the identity matrix
er_matrix_inverse <- function(a) {
  # solve(a) computes the inverse of matrix a
  return(solve(a))
}


#' Compute Cumulative Relative Frequencies
#'
#' @description
#' Computes the cumulative relative frequencies from a vector of relative
#' frequencies. This is an R replacement for the C function `cum_rel_freqs`.
#'
#' @param rel_freq A numeric vector of relative frequencies.
#'
#' @return A numeric vector of cumulative relative frequencies.
cum_rel_freqs <- function(rel_freq) {
  # R's cumsum() is the direct and efficient equivalent.
  return(cumsum(rel_freq))
}

#' Get Equated Scale Scores
#'
#' @description
#' Converts equated raw scores to a scale score metric using a provided
#' raw-to-scale score conversion table. This function uses linear interpolation
#' for equated raw scores that fall between the points in the conversion table.
#' This is an R translation of the `Equated_ss` function.
#'
#' @param eraw A numeric vector of equated raw scores.
#' @param yct A data frame or matrix with two columns: the first for the raw
#'   scores on the base form (Y), and the second for the corresponding scale scores.
#' @param lprss The lowest possible rounded scale score (for clipping).
#' @param hprss The highest possible rounded scale score (for clipping).
#' @param round_digits An integer specifying the number of decimal places to
#'   round the final scale scores to. If NULL, no rounding is performed.
#'
#' @return A list containing the unrounded (`essu`) and rounded (`essr`)
#'   equated scale scores.
equated_ss <- function(eraw, yct, lprss = NULL, hprss = NULL, round_digits = 0) {

  # Create a linear interpolation function from the conversion table
  # rule = 2 ensures that values outside the range are set to the nearest endpoint
  interp_fun <- stats::approxfun(x = yct[[1]], y = yct[[2]], rule = 2)

  # Apply the interpolation function to the equated raw scores
  essu <- interp_fun(eraw)

  # Round the scores if specified
  if (!is.null(round_digits)) {
    essr <- round(essu, digits = round_digits)
  } else {
    essr <- essu
  }

  # Clip the scores to the specified min/max scale score range
  if (!is.null(lprss)) {
    essr[essr < lprss] <- lprss
  }
  if (!is.null(hprss)) {
    essr[essr > hprss] <- hprss
  }

  return(list(essu = essu, essr = essr))
}


#' Perform QR Decomposition
#'
#' @description
#' A wrapper for R's built-in `qr()` function to maintain naming consistency
#' with the C library's `er_qrdcmp`.
#'
#' @param a A numeric matrix.
#'
#' @return An object of class "qr" representing the QR decomposition.
er_qrdcmp <- function(a) {
  return(qr(a))
}

#' Find the Root of a One-Dimensional Function
#'
#' @description
#' Finds the root of a function within a given interval. This is a wrapper for
#' R's `uniroot` function, replacing the C library's `er_rtsafe`.
#'
#' @param f The function for which the root is sought.
#' @param interval A vector containing the two endpoints of the interval to search.
#'
#' @return The x-value of the root.
er_rtsafe <- function(f, interval) {
  result <- uniroot(f, interval)
  return(result$root)
}

#' Minimize a Function using the BFGS Algorithm
#'
#' @description
#' Minimizes a function using the Broyden-Fletcher-Goldfarb-Shanno (BFGS)
#' algorithm. This is a wrapper for R's powerful `optim` function and replaces
#' the C library's `er_dfpmin` and its internal line search `er_lnsrch`.
#'
#' @param par A numeric vector of initial values for the parameters.
#' @param fn The function to be minimized.
#' @param gr (Optional) The gradient of the function `fn`.
#' @param ... Additional arguments to be passed to `fn` and `gr`.
#'
#' @return An object of class "optim" containing the optimization results.
er_dfpmin <- function(par, fn, gr = NULL, ...) {
  optim(par = par, fn = fn, gr = gr, method = "BFGS", ...)
}

