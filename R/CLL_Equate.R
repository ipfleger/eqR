# Title: R Functions for Continuized Log-Linear (CLL) Equating (Fully Vectorized)
#
# Description:
# This script contains fully vectorized R translations of the functions from
# 'CLL_Equate.c'. All high-level equating functions have been optimized to
# use matrix operations instead of loops for significant performance gains.

#' Vectorized Bivariate Exponential Polynomial Calculation
#'
#' @description
#' The core computational engine for the log-linear model's density function.
#' This version is highly optimized for speed using matrix operations.
biv_exp_polynomial_vec <- function(bivar, x, y) {
  full_matrix <- matrix(0, nrow = length(x), ncol = length(bivar$Beta))
  if (bivar$cu > 0) for (i in seq_along(bivar$cu_indices)) full_matrix[, bivar$cu_indices[i]] <- x^i
  if (bivar$cv > 0) for (i in seq_along(bivar$cv_indices)) full_matrix[, bivar$cv_indices[i]] <- y^i
  if (bivar$cuv > 0) for (i in seq_along(bivar$cuv_indices)) full_matrix[, bivar$cuv_indices[i]] <- x^bivar$cpm[i, 1] * y^bivar$cpm[i, 2]
  exponent <- full_matrix %*% bivar$Beta
  return(as.vector(exp(bivar$ap + exponent)))
}

#' Gaussian Quadrature for Numerical Integration (Fixed-Point)
#'
#' @description
#' Approximates a definite integral using a 20-point Legendre-Gauss fixed-point rule.
gaussian_quadrature_vec <- function(f, lower, upper, ...) {
  nodes <- c(-0.9931285991850949, -0.9639719272779138, -0.9122344282513260, -0.8391169718222188, -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154196, -0.2277858511416451, -0.0765265211334973, 0.0765265211334973, 0.2277858511416451, 0.3737060887154196, 0.5108670019508271, 0.6360536807265150, 0.7463319064601508, 0.8391169718222188, 0.9122344282513260, 0.9639719272779138, 0.9931285991850949)
  weights <- c(0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767047, 0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183821, 0.1491729864726037, 0.1527533871307258, 0.1527533871307258, 0.1491729864726037, 0.1420961093183821, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404, 0.0832767415767047, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521)
  t_nodes <- 0.5 * ((upper - lower) * nodes + (upper + lower))
  f_values <- f(t_nodes)
  integral <- 0.5 * (upper - lower) * sum(weights * f_values)
  attr(integral, "nodes") <- t_nodes
  attr(integral, "weights") <- weights
  return(integral)
}

#' Fully Vectorized Integrand Function
#'
#' @description
#' A helper function that computes the integral over the y-dimension for all
#' x-values at once, avoiding slow R loops.
integrand_fully_vectorized <- function(x_vec, bivar, miny, maxy) {
  quad_info <- gaussian_quadrature_vec(function(y) y, lower = miny, upper = maxy)
  y_nodes <- attr(quad_info, "nodes")
  y_weights <- attr(quad_info, "weights")
  grid <- expand.grid(x = x_vec, y = y_nodes)
  poly_values <- biv_exp_polynomial_vec(bivar, grid$x, grid$y)
  value_matrix <- matrix(poly_values, nrow = length(x_vec), ncol = length(y_nodes), byrow = FALSE)
  integral_values <- 0.5 * (maxy - miny) * (value_matrix %*% y_weights)
  return(as.vector(integral_values))
}

#' Calculate Marginal CDF for X
cll_marginal_cdf_x <- function(x, bivar, nc) {
  minx <- bivar$minx - 0.5
  miny <- bivar$minv - 0.5
  maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5
  integrand <- function(x_vec) integrand_fully_vectorized(x_vec, bivar, miny, maxy)
  numerator <- gaussian_quadrature_vec(integrand, lower = minx, upper = x)
  return(numerator / nc)
}

#' Calculate Marginal CDF for Y
cll_marginal_cdf_y <- function(y, bivar, nc) {
  minx <- bivar$minx - 0.5
  maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  miny <- bivar$minv - 0.5
  integrand <- function(y_vec) {
    # This still needs to be an inner integration
    sapply(y_vec, function(yi) {
      integrand_x <- function(x_vec) biv_exp_polynomial_vec(bivar, x_vec, rep(yi, length(x_vec)))
      gaussian_quadrature_vec(integrand_x, lower = minx, upper = maxx)
    })
  }
  numerator <- gaussian_quadrature_vec(integrand, lower = miny, upper = y)
  return(numerator / nc)
}

#' Calculate Inverse Marginal CDF for X
cll_marginal_inverse_cdf_x <- function(p, bivar, nc) {
  minx <- bivar$minx - 0.5
  maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  f_to_solve <- function(x) cll_marginal_cdf_x(x, bivar, nc) - p
  uniroot(f_to_solve, interval = c(minx, maxx))$root
}

#' Calculate Inverse Marginal CDF for Y
cll_marginal_inverse_cdf_y <- function(p, bivar, nc) {
  miny <- bivar$minv - 0.5
  maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5
  f_to_solve <- function(y) cll_marginal_cdf_y(y, bivar, nc) - p
  uniroot(f_to_solve, interval = c(miny, maxy))$root
}

# --- High-Level Equating Functions (All Vectorized) ---

cll_equate_sg <- function(bivar) {
  scoresx <- bivar$minx + (0:(bivar$nsx - 1)) * bivar$incx
  minx <- bivar$minx - 0.5; maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  miny <- bivar$minv - 0.5; maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5

  integrand_for_nc <- function(x_vec) integrand_fully_vectorized(x_vec, bivar, miny, maxy)
  nc <- gaussian_quadrature_vec(integrand_for_nc, lower = minx, upper = maxx)

  cdfx <- sapply(scoresx, function(x_val) cll_marginal_cdf_x(x_val, bivar, nc))
  sapply(cdfx, function(p) cll_marginal_inverse_cdf_y(p, bivar, nc))
}

cll_equate_cb <- function(bivar1, bivar2, wtsx, wtsy) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  integrand1 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar1, miny, maxy)
  integrand2 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar2, miny, maxy)
  nc1 <- gaussian_quadrature_vec(integrand1, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature_vec(integrand2, lower = minx, upper = maxx)

  cdfx <- sapply(scoresx, function(x_val) wtsx * cll_marginal_cdf_x(x_val, bivar1, nc1) + (1 - wtsx) * cll_marginal_cdf_x(x_val, bivar2, nc2))

  cdf_y_combined <- function(y) {
    wtsy * cll_marginal_cdf_y(y, bivar1, nc1) + (1 - wtsy) * cll_marginal_cdf_y(y, bivar2, nc2)
  }

  inverse_cdf_y_combined <- function(p) {
    uniroot(function(y) cdf_y_combined(y) - p, interval = c(miny, maxy))$root
  }

  sapply(cdfx, inverse_cdf_y_combined)
}

cll_equate_neat_ps <- function(bivar1, bivar2, wts) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  integrand1 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar1, miny, maxy)
  integrand2 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar2, miny, maxy)
  nc1 <- gaussian_quadrature_vec(integrand1, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature_vec(integrand2, lower = minx, upper = maxx)

  cdf_x_ps <- function(x) wts * cll_marginal_cdf_x(x, bivar1, nc1) + (1 - wts) * cll_marginal_cdf_x(x, bivar2, nc2)
  cdf_y_ps <- function(y) wts * cll_marginal_cdf_y(y, bivar1, nc1) + (1 - wts) * cll_marginal_cdf_y(y, bivar2, nc2)

  cdfx <- sapply(scoresx, cdf_x_ps)

  inverse_cdf_y_ps <- function(p) {
    uniroot(function(y) cdf_y_ps(y) - p, interval = c(miny, maxy))$root
  }

  sapply(cdfx, inverse_cdf_y_ps)
}

cll_equate_neat_chn <- function(bivar1, bivar2) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  integrand1 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar1, miny, maxy)
  integrand2 <- function(x_vec) integrand_fully_vectorized(x_vec, bivar2, miny, maxy)
  nc1 <- gaussian_quadrature_vec(integrand1, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature_vec(integrand2, lower = minx, upper = maxx)

  cdfx1 <- sapply(scoresx, function(x) cll_marginal_cdf_x(x, bivar = bivar1, nc1))
  v_scores <- sapply(cdfx1, function(p) cll_marginal_inverse_cdf_y(p, bivar = bivar1, nc1))
  cdfv2 <- sapply(v_scores, function(v) cll_marginal_cdf_y(v, bivar2, nc2))
  equated_scores <- sapply(cdfv2, function(p) cll_marginal_inverse_cdf_x(p, bivar2, nc2))

  return(equated_scores)
}

#' High-Level Wrapper for all CLL Equating Designs
equate_cll <- function(design, ...) {
  args <- list(...)
  switch(design,
         "SG" = cll_equate_sg_cpp(args$bivar),
         "CB" = cll_equate_cb(args$bivar1, args$bivar2, args$wtsx, args$wtsy),
         "NEAT_PS" = cll_equate_neat_ps(args$bivar1, args$bivar2, args$wts),
         "NEAT_CHN" = cll_equate_neat_chn_rcpp(args$bivar1, args$bivar2),
         stop("Invalid design specified.")
  )
}
# ------------------------------------------------------------------------------
# Above this point was the first draft
#' Perform CLL Equating (NEAT Chained) using Rcpp for Speed
#'
#' @description
#' A high-performance version of `cll_equate_neat_chn` that calls the compiled
#' C++ functions for the most intensive calculations.
#'
#' @param bivar1, bivar2 The `bivar` objects for the two forms.
#'
#' @return A numeric vector of equated scores.
cll_equate_neat_chn_rcpp <- function(bivar1, bivar2) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  # Pre-calculate Gaussian quadrature info and add it to bivar objects
  # This avoids recalculating it inside the C++ loop
  quad_info_y <- gaussian_quadrature_vec(function(y) y, lower = miny, upper = maxy)
  bivar1$quad_info <- list(nodes = attr(quad_info_y, "nodes"), weights = attr(quad_info_y, "weights"))
  bivar2$quad_info <- list(nodes = attr(quad_info_y, "nodes"), weights = attr(quad_info_y, "weights"))

  # Define the integrand function that calls our C++ code
  integrand1_cpp <- function(x_vec) integrand_fully_vectorized_cpp(bivar1, x_vec, miny, maxy)
  integrand2_cpp <- function(x_vec) integrand_fully_vectorized_cpp(bivar2, x_vec, miny, maxy)

  # Calculate normalization constants using the C++ functions
  nc1 <- gaussian_quadrature_vec(integrand1_cpp, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature_vec(integrand2_cpp, lower = minx, upper = maxx)

  # The rest of the logic remains in R, as these are not the primary bottlenecks
  cdfx1 <- sapply(scoresx, function(x) cll_marginal_cdf_x(x, bivar = bivar1, nc1))
  v_scores <- sapply(cdfx1, function(p) cll_marginal_inverse_cdf_y(p, bivar = bivar1, nc1))
  cdfv2 <- sapply(v_scores, function(v) cll_marginal_cdf_y(v, bivar2, nc2))
  equated_scores <- sapply(cdfv2, function(p) cll_marginal_inverse_cdf_x(p, bivar2, nc2))

  return(equated_scores)
}

# --------------------------------------------------------------------
# CORRECTED VERSION - 2025-07-24 (v2)

# Helper function to get Gauss-Legendre quadrature nodes and weights
.get_gauss_legendre_params <- function() {
  # 20-point Legendre-Gauss nodes and weights for the interval [-1, 1]
  nodes <- c(
    -0.9931285991850949, -0.9639719272779138, -0.9122344282513260,
    -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
    -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
    -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
    0.3737060887154196,  0.5108670019508271,  0.6360536807265150,
    0.7463319064601508,  0.8391169718222188,  0.9122344282513260,
    0.9639719272779138,  0.9931285991850949
  )
  weights <- c(
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
    0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
    0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
    0.1527533871307259, 0.1527533871307259, 0.1491729864726037,
    0.1420961093183821, 0.1316886384491766, 0.1181945319615184,
    0.1019301198172404, 0.0832767415767048, 0.0626720483341091,
    0.0406014298003869, 0.0176140071391521
  )
  return(list(nodes = nodes, weights = weights))
}


#' Numerically Integrate Using C++ Backend (Vectorized)
#'
#' @description
#' This function performs a double integral using Gauss-Legendre quadrature.
#' It calls the C++ function `integrand_fully_vectorized_cpp` to evaluate
#' the inner function efficiently over a grid.
#'
#' @param bivar The bivariate model object.
#' @param lower_x Lower integration limit for X.
#' @param upper_x Upper integration limit for X.
#' @param lower_y Lower integration limit for Y.
#' @param upper_y Upper integration limit for Y.
#' @return The value of the double integral.
gaussian_quadrature_cpp <- function(bivar, lower_x, upper_x, lower_y, upper_y) {
  gl <- .get_gauss_legendre_params()

  # Transform nodes for both X and Y dimensions
  x_nodes <- 0.5 * (upper_x - lower_x) * gl$nodes + 0.5 * (upper_x + lower_x)
  y_nodes <- 0.5 * (upper_y - lower_y) * gl$nodes + 0.5 * (upper_y + lower_y)

  # --- START: NEW CODE to fix the crash ---
  # 1. Generate x_powers and y_powers vectors based on the bivar object
  x_powers <- c(0, 1:bivar$cu, rep(0, bivar$cv), bivar$cpm[, 1])
  y_powers <- c(0, rep(0, bivar$cu), 1:bivar$cv, bivar$cpm[, 2])

  # 2. Extract the model coefficients from bivar$Beta
  model_params <- as.vector(bivar$Beta)
  # --- END: NEW CODE ---

  # Call the C++ function with the CORRECT parameters
  integrand_values <- integrand_fully_vectorized_cpp(
    x = x_nodes,
    y_nodes = y_nodes,
    bivar_params = model_params, # Use the corrected parameter
    x_powers = x_powers,         # Use the generated powers
    y_powers = y_powers          # Use the generated powers
  )

  # Perform the integration using matrix multiplication
  inner_integral <- t(integrand_values) %*% (gl$weights * 0.5 * (upper_y - lower_y))
  outer_integral <- sum(inner_integral * gl$weights * 0.5 * (upper_x - lower_x))

  return(outer_integral)
}

#' Calculate Marginal CDF for X (C++ Accelerated)
#'
#' @description
#' Calculates the cumulative distribution function P(X <= x_val).
cll_marginal_cdf_x_cpp <- function(x_val, bivar, nc) {
  minx <- bivar$minx - 0.5
  miny <- bivar$minv - 0.5
  maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5

  # Integrate from minx to x_val
  numerator <- gaussian_quadrature_cpp(bivar, minx, x_val, miny, maxy)
  return(numerator / nc)
}

#' Calculate Marginal CDF for Y (C++ Accelerated)
#'
#' @description
#' Calculates the cumulative distribution function P(Y <= y_val).
cll_marginal_cdf_y_cpp <- function(y_val, bivar, nc) {
  minx <- bivar$minx - 0.5
  maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  miny <- bivar$minv - 0.5

  # Integrate from miny to y_val
  numerator <- gaussian_quadrature_cpp(bivar, minx, maxx, miny, y_val)
  return(numerator / nc)
}


#' Calculate Inverse Marginal CDF for Y (C++ Accelerated)
#'
#' @description
#' Finds the score y such that P(Y <= y) = p.
cll_marginal_inverse_cdf_y_cpp <- function(p, bivar, nc) {
  miny <- bivar$minv - 0.5
  maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5

  objective_fn <- function(y_val) {
    cll_marginal_cdf_y_cpp(y_val, bivar, nc) - p
  }

  result <- try(stats::uniroot(objective_fn, interval = c(miny, maxy))$root, silent = TRUE)
  if (inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}

#' CLL Single-Group Equating (C++ Accelerated)
#'
#' @description
#' This is the main wrapper function for single-group CLL equating.
#'
#' @param bivar A bivariate model object created by `smooth_bll`.
#' @return A vector of equated Y scores corresponding to the input X scores.
#' @export
cll_equate_sg_cpp <- function(bivar) {
  scoresx <- bivar$minx + (0:(bivar$nsx - 1)) * bivar$incx
  minx <- bivar$minx - 0.5
  maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  miny <- bivar$minv - 0.5
  maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5

  nc <- gaussian_quadrature_cpp(bivar = bivar, lower_x = minx, upper_x = maxx, lower_y = miny, upper_y = maxy)
  if (is.na(nc) || nc == 0) {
    warning("Normalizing constant is zero or NA. Cannot proceed with equating.")
    return(rep(NA, bivar$nsx))
  }

  cdfx <- sapply(scoresx, function(x_val) cll_marginal_cdf_x_cpp(x_val, bivar, nc))

  sapply(cdfx, function(p) {
    if (is.na(p)) return(NA)
    cll_marginal_inverse_cdf_y_cpp(p, bivar, nc)
  })
}
