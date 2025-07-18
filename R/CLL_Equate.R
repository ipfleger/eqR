# Title: R Functions for Continuized Log-Linear (CLL) Equating
#
# Description:
# This script contains R translations of the functions from 'CLL_Equate.c'.
# These functions implement the continuized log-linear equating method
# described by Holland and Thayer (2004) and referenced in von Davier et al. (2004).

# Source the file containing previously converted utility functions.
# Note: You may need to change this file path.
# source("equating_recipes_utils.R")
# source("r_loglinear_utils.R") # For smoothing functions

#' Evaluate the Exponential of a Polynomial (Vectorized)
#'
#' @description
#' The core component of the log-linear model's density function.
#' Translates `ExpPolynomial`. This version is vectorized to work with `integrate`.
exp_polynomial <- function(x_vec, params) {
  sapply(x_vec, function(x) {
    powers <- 0:(length(params) - 1)
    s <- sum(params * (x^powers))
    return(exp(s))
  })
}

#' Numerical Integration using R's `integrate`
#'
#' @description
#' A general-purpose numerical integration function to replace the specific
#' `GaussianQuadrature` functions in the C code.
gaussian_quadrature <- function(f, lower, upper, ...) {
  result <- try(stats::integrate(f, lower, upper, ...)$value, silent = TRUE)
  if (inherits(result, "try-error")) {
    warning("Integration failed, returning 0.")
    return(0)
  }
  return(result)
}

  #' Gaussian Quadrature for Numerical Integration (Final Version)
  #'
  #' @description
  #' Approximates a definite integral using a 20-point Legendre-Gauss fixed-point
  #' rule. This function is a drop-in replacement for the previous `stats::integrate`
  #' version, designed for speed and vectorization.
  #'
  #' @param f The function to integrate.
  #' @param lower The lower limit of integration.
  #' @param upper The upper limit of integration.
  #' @param ... Catches any additional arguments from the calling function to
  #'   ensure signature compatibility. These are not used by this fixed-point method.
  #'
  #' @return The approximate value of the integral as a single numeric value.
  gaussian_quadrature_vec <- function(f, lower, upper, ...) {
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
      0.0832767415767047, 0.1019301198172404, 0.1181945319615184,
      0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
      0.1527533871307258, 0.1527533871307258, 0.1491729864726037,
      0.1420961093183821, 0.1316886384491766, 0.1181945319615184,
      0.1019301198172404, 0.0832767415767047, 0.0626720483341091,
      0.0406014298003869, 0.0176140071391521
    )

    # Transform nodes to the [lower, upper] interval
    t_nodes <- 0.5 * ((upper - lower) * nodes + (upper + lower))

    f_values <- f(t_nodes)

    integral <- 0.5 * (upper - lower) * sum(weights * f_values)

    # Attach nodes and weights for internal use by the vectorized integrand
    attr(integral, "nodes") <- t_nodes
    attr(integral, "weights") <- weights

    # Return the numeric value directly
    return(integral)
  }
#' Continuized Log-Linear PDF for Equivalent Groups
#'
#' @description Translates `CLLEGPdf`.
cll_eg_pdf <- function(x, min, max, params) {
  nc <- gaussian_quadrature(exp_polynomial, lower = min, upper = max, params = params)
  if (nc == 0) return(0)
  return(exp_polynomial(x, params) / nc)
}

#' Continuized Log-Linear CDF for Equivalent Groups
#'
#' @description Translates `CLLEGCdf`.
cll_eg_cdf <- function(x, min, max, params) {
  nc <- gaussian_quadrature(exp_polynomial, lower = min, upper = max, params = params)
  if (nc == 0) return(ifelse(x >= max, 1, 0))
  cdf_val <- gaussian_quadrature(exp_polynomial, lower = min, upper = x, params = params)
  return(cdf_val / nc)
}

#' Inverse of the Continuized Log-Linear CDF
#'
#' @description Translates `CLLInverseCdf`.
cll_inverse_cdf <- function(p, min, max, params) {
  f_to_solve <- function(x) cll_eg_cdf(x, min, max, params) - p
  if (f_to_solve(min) * f_to_solve(max) > 1e-8) { # Allow for floating point inaccuracies
    if(abs(f_to_solve(min)) < 1e-6) return(min)
    if(abs(f_to_solve(max)) < 1e-6) return(max)
    warning("Root not bracketed in cll_inverse_cdf. Check inputs.")
    return(NA)
  }
  uniroot(f_to_solve, interval = c(min, max))$root
}

#' Bivariate Exponential Polynomial
#'
#' @description Translates `BivExpPolynomial`. This version is NOT vectorized.
biv_exp_polynomial <- function(bivar_list, x, y) {
  # Unpack parameters from the list
  nx <- bivar_list$cu
  ny <- bivar_list$cv
  nxbyy <- bivar_list$cuv

  s <- bivar_list$ap
  if (nx > 0) s <- s + sum(bivar_list$Beta[1:nx] * (x^(1:nx)))
  if (ny > 0) s <- s + sum(bivar_list$Beta[(nx + 1):(nx + ny)] * (y^(1:ny)))
  if (nxbyy > 0) {
    for (j in 1:nxbyy) {
      s <- s + bivar_list$Beta[nx + ny + j] * (x^bivar_list$cpm[j, 1]) * (y^bivar_list$cpm[j, 2])
    }
  }
  return(exp(s))
}

#' Calculate Marginal PDF for X from a Bivariate Distribution
cll_marginal_pdf_x <- function(x_vec, bivar, nc) {
  sapply(x_vec, function(x) {
    miny <- bivar$minv - 0.5
    maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5
    integrand <- function(y_vec_inner) {
      sapply(y_vec_inner, function(y) biv_exp_polynomial(bivar, x, y))
    }
    val <- gaussian_quadrature(integrand, lower = miny, upper = maxy)
    return(val / nc)
  })
}

#' Calculate Marginal PDF for Y from a Bivariate Distribution
cll_marginal_pdf_y <- function(y_vec, bivar, nc) {
  sapply(y_vec, function(y) {
    minx <- bivar$minx - 0.5
    maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
    integrand <- function(x_vec_inner) {
      sapply(x_vec_inner, function(x) biv_exp_polynomial(bivar, x, y))
    }
    val <- gaussian_quadrature(integrand, lower = minx, upper = maxx)
    return(val / nc)
  })
}

#' Calculate Marginal CDF for X
cll_marginal_cdf_x <- function(x, bivar, nc) {
  minx <- bivar$minx - 0.5
  gaussian_quadrature(cll_marginal_pdf_x, lower = minx, upper = x, bivar = bivar, nc = nc)
}

#' Calculate Marginal CDF for Y
cll_marginal_cdf_y <- function(y, bivar, nc) {
  miny <- bivar$minv - 0.5
  gaussian_quadrature(cll_marginal_pdf_y, lower = miny, upper = y, bivar = bivar, nc = nc)
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

#' Perform CLL Equating for the Single Group Design
cll_equate_sg <- function(bivar) {
  scoresx <- bivar$minx + (0:(bivar$nsx - 1)) * bivar$incx
  minx <- bivar$minx - 0.5; maxx <- bivar$minx + (bivar$nsx - 1) * bivar$incx + 0.5
  miny <- bivar$minv - 0.5; maxy <- bivar$minv + (bivar$nsv - 1) * bivar$incv + 0.5

  integrand_for_nc <- function(x_vec) {
    sapply(x_vec, function(xi) {
      integrand_y <- function(y_vec) sapply(y_vec, function(yi) biv_exp_polynomial(bivar, xi, yi))
      gaussian_quadrature(integrand_y, lower = miny, upper = maxy)
    })
  }
  nc <- gaussian_quadrature(integrand_for_nc, lower = minx, upper = maxx)

  cdfx <- sapply(scoresx, function(x_val) cll_marginal_cdf_x(x_val, bivar, nc))
  sapply(cdfx, function(p) cll_marginal_inverse_cdf_y(p, bivar, nc))
}

#' Perform CLL Equating for the Counter-Balanced Design
cll_equate_cb <- function(bivar1, bivar2, wtsx, wtsy) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  integrand_for_nc1 <- function(x_vec) sapply(x_vec, function(xi) gaussian_quadrature(function(yi_vec) sapply(yi_vec, function(yi) biv_exp_polynomial(bivar1, xi, yi)), lower = miny, upper = maxy))
  integrand_for_nc2 <- function(x_vec) sapply(x_vec, function(xi) gaussian_quadrature(function(yi_vec) sapply(yi_vec, function(yi) biv_exp_polynomial(bivar2, xi, yi)), lower = miny, upper = maxy))
  nc1 <- gaussian_quadrature(integrand_for_nc1, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature(integrand_for_nc2, lower = minx, upper = maxx)

  cdfx <- sapply(scoresx, function(x_val) wtsx * cll_marginal_cdf_x(x_val, bivar1, nc1) + (1 - wtsx) * cll_marginal_cdf_x(x_val, bivar2, nc2))

  cll_marginal_cdf_cb_y <- function(y_vec) wtsy * cll_marginal_pdf_y(y_vec, bivar1, nc1) + (1 - wtsy) * cll_marginal_pdf_y(y_vec, bivar2, nc2)
  cll_marginal_inverse_cdf_cb_y <- function(p) {
    uniroot(function(y) gaussian_quadrature(cll_marginal_cdf_cb_y, lower=miny, upper=y) - p, interval = c(miny, maxy))$root
  }

  sapply(cdfx, cll_marginal_inverse_cdf_cb_y)
}

#' Perform CLL Equating for NEAT Post-Stratification Design
cll_equate_neat_ps <- function(bivar1, bivar2, wts) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  integrand_for_nc1 <- function(x_vec) sapply(x_vec, function(xi) gaussian_quadrature(function(yi_vec) sapply(yi_vec, function(yi) biv_exp_polynomial(bivar1, xi, yi)), lower = miny, upper = maxy))
  integrand_for_nc2 <- function(x_vec) sapply(x_vec, function(xi) gaussian_quadrature(function(yi_vec) sapply(yi_vec, function(yi) biv_exp_polynomial(bivar2, xi, yi)), lower = miny, upper = maxy))
  nc1 <- gaussian_quadrature(integrand_for_nc1, lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature(integrand_for_nc2, lower = minx, upper = maxx)

  cll_neat_ps_marginal_pdf <- function(score_vec, is_x) {
    pdf1 <- if(is_x) cll_marginal_pdf_x(score_vec, bivar1, nc1) else cll_marginal_pdf_y(score_vec, bivar1, nc1)
    pdf2 <- if(is_x) cll_marginal_pdf_x(score_vec, bivar2, nc2) else cll_marginal_pdf_y(score_vec, bivar2, nc2)
    return(wts * pdf1 + (1 - wts) * pdf2)
  }

  cll_neat_ps_marginal_cdf <- function(score, is_x) {
    lower_bound <- if(is_x) minx else miny
    gaussian_quadrature(function(s) cll_neat_ps_marginal_pdf(s, is_x), lower = lower_bound, upper = score)
  }

  cdfx <- sapply(scoresx, function(x) cll_neat_ps_marginal_cdf(x, is_x = TRUE))

  cll_neat_ps_inverse_cdf <- function(p) {
    uniroot(function(y) cll_neat_ps_marginal_cdf(y, is_x=FALSE) - p, interval = c(miny, maxy))$root
  }

  sapply(cdfx, cll_neat_ps_inverse_cdf)
}

#' Vectorized Bivariate Exponential Polynomial Calculation (Corrected)
#'
#' @description
#' A corrected, vectorized version that uses pre-calculated indices from the
#' bivar object to correctly and efficiently build the design matrix.
#'
#' @param bivar A `bivar` object from the smoothing process.
#' @param x A vector of scores for the first variable (u).
#' @param y A vector of scores for the second variable (v).
#'
#' @return A numeric vector of results from the exponential polynomial calculation.
biv_exp_polynomial_vec <- function(bivar, x, y) {
  # Initialize the full design matrix
  full_matrix <- matrix(0, nrow = length(x), ncol = length(bivar$Beta))

  # Populate the matrix using the pre-calculated indices
  # This is both faster and less error-prone.
  if (bivar$cu > 0) {
    for (i in seq_along(bivar$cu_indices)) {
      full_matrix[, bivar$cu_indices[i]] <- x^i
    }
  }
  if (bivar$cv > 0) {
    for (i in seq_along(bivar$cv_indices)) {
      full_matrix[, bivar$cv_indices[i]] <- y^i
    }
  }
  if (bivar$cuv > 0) {
    for (i in seq_along(bivar$cuv_indices)) {
      full_matrix[, bivar$cuv_indices[i]] <- x^bivar$cpm[i, 1] * y^bivar$cpm[i, 2]
    }
  }

  # Calculate the exponent with a single matrix multiplication
  exponent <- full_matrix %*% bivar$Beta

  return(as.vector(exp(bivar$ap + exponent)))
}


#' Perform CLL Equating for NEAT Chained Design
cll_equate_neat_chn <- function(bivar1, bivar2) {
  scoresx <- bivar1$minx + (0:(bivar1$nsx - 1)) * bivar1$incx
  minx <- bivar1$minx - 0.5; maxx <- bivar1$minx + (bivar1$nsx - 1) * bivar1$incx + 0.5
  miny <- bivar1$minv - 0.5; maxy <- bivar1$minv + (bivar1$nsv - 1) * bivar1$incv + 0.5

  # This is the new, fully vectorized function that computes the integrals for ALL x values at once.
  integrand_fully_vectorized <- function(x_vec, bivar) {
    # Get the integration nodes and weights for the y-dimension.
    # We call gaussian_quadrature on a dummy function just to get the attributes.
    quad_info <- gaussian_quadrature_vec(function(y) y, lower = miny, upper = maxy)
    y_nodes <- attr(quad_info, "nodes")
    y_weights <- attr(quad_info, "weights")

    # 1. Create a grid of all (x, y) combinations needed for the integration.
    grid <- expand.grid(x = x_vec, y = y_nodes)

    # 2. Call the expensive polynomial function only ONCE on the entire grid.
    poly_values <- biv_exp_polynomial_vec(bivar, grid$x, grid$y)

    # 3. Reshape the results into a matrix where rows correspond to x and columns to y.
    value_matrix <- matrix(poly_values, nrow = length(x_vec), ncol = length(y_nodes), byrow = FALSE)

    # 4. Perform the integration for all x values simultaneously using fast matrix multiplication.
    integral_values <- 0.5 * (maxy - miny) * (value_matrix %*% y_weights)

    return(as.vector(integral_values))
  }

  # The outer call to gaussian_quadrature now uses our fast, vectorized integrand.
  nc1 <- gaussian_quadrature_vec(function(x_vec) integrand_fully_vectorized(x_vec, bivar1), lower = minx, upper = maxx)
  nc2 <- gaussian_quadrature_vec(function(x_vec) integrand_fully_vectorized(x_vec, bivar2), lower = minx, upper = maxx)

  # --- The rest of the function remains the same ---
  cdfx1 <- sapply(scoresx, function(x) cll_marginal_cdf_x(x, bivar = bivar1, nc1))
  v_scores <- sapply(cdfx1, function(p) cll_marginal_inverse_cdf_y(p, bivar = bivar1, nc1))
  cdfv2 <- sapply(v_scores, function(v) cll_marginal_cdf_y(v, bivar2, nc2))
  equated_scores <- sapply(cdfv2, function(p) cll_marginal_inverse_cdf_x(p, bivar2, nc2))

  return(equated_scores)
}
#' Perform High-Level Continuized Log-Linear Equating
#'
#' @description
#' This function is a high-level wrapper that performs continuized log-linear
#' equating for various designs by calling the appropriate helper functions.
#' It replaces the C functions `Wrapper_SC`, `Wrapper_CC`, and `Wrapper_BC`.
#'
#' @param design A string specifying the equating design. One of "SG" (Single
#'   Group), "CB" (Counter-Balanced), "NEAT_PS" (NEAT Post-Stratification), or
#'   "NEAT_CHN" (NEAT Chained).
#' @param ... Additional arguments required by the specific design, such as
#'   `bivar1`, `bivar2`, `wts`, etc.
#'
#' @return A vector of equated scores.
equate_cll <- function(design, ...) {# w1 = 1, # We may want to expose this one specifically. Three common values (according to Gemini) are 1, .5, or -1. I'm not entirely sure what they do

  args <- list(...)

  switch(design,
         "SG" = {
           if (!"bivar" %in% names(args)) stop("Single Group design requires 'bivar' argument.")
           return(cll_equate_sg(args$bivar))
         },
         "CB" = {
           if (!all(c("bivar1", "bivar2", "wtsx", "wtsy") %in% names(args))) {
             stop("Counter-Balanced design requires 'bivar1', 'bivar2', 'wtsx', and 'wtsy'.")
           }
           return(cll_equate_cb(args$bivar1, args$bivar2, args$wtsx, args$wtsy))
         },
         "NEAT_PS" = {
           if (!all(c("bivar1", "bivar2", "wts") %in% names(args))) {
             stop("NEAT Post-Stratification design requires 'bivar1', 'bivar2', and 'wts'.")
           }
           return(cll_equate_neat_ps(args$bivar1, args$bivar2, args$wts))
         },
         "NEAT_CHN" = {
           if (!all(c("bivar1", "bivar2") %in% names(args))) {
             stop("NEAT Chained design requires 'bivar1' and 'bivar2'.")
           }
           return(cll_equate_neat_chn(args$bivar1, args$bivar2))
         },
         stop("Invalid design specified.")
  )
}
