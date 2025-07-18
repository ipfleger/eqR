# Title: R Functions for Kernel Equating and Standard Error Estimation
#
# Description:
# This script contains R translations of the kernel equating functions from
# the "Equating Recipes" C library, found in 'kernel_Equate.c'. These
# functions are used to perform kernel smoothing of score distributions,
# find the optimal bandwidth parameter, perform the equating, and calculate
# analytic standard errors using the delta method.

# --- Core Kernel Density and Distribution Functions ---

#' Kernel Smoothed Continuous Probability Density Function (PDF)
#'
#' @description
#' Computes the kernel smoothed continuous PDF for a given score `x`.
#'
#' @param x A numeric value at which to evaluate the PDF.
#' @param scores A numeric vector of the discrete score points.
#' @param rel_freq A numeric vector of the relative frequencies for each score.
#' @param hx A numeric value for the bandwidth parameter.
#' @return The value of the kernel smoothed PDF at `x`.
kernel_continu_pdf <- function(x, scores, rel_freq, hx) {
  mu <- sum(scores * rel_freq)
  sigma_sq <- sum(scores^2 * rel_freq) - mu^2
  if (sigma_sq <= 0) return(0)
  ax <- sqrt(sigma_sq / (sigma_sq + hx^2))
  z_scores <- (x - ax * scores - (1 - ax) * mu) / (ax * hx)
  pdf_val <- sum(rel_freq * dnorm(z_scores) / (ax * hx))
  return(pdf_val)
}

#' Kernel Smoothed Continuous Cumulative Distribution Function (CDF)
#'
#' @description
#' Computes the kernel smoothed continuous CDF for a given score `x`.
#'
#' @param x A numeric value at which to evaluate the CDF.
#' @param scores A numeric vector of the discrete score points.
#' @param rel_freq A numeric vector of the relative frequencies for each score.
#' @param hx A numeric value for the bandwidth parameter.
#' @return The value of the kernel smoothed CDF at `x`.
kernel_continu_cdf <- function(x, scores, rel_freq, hx) {
  mu <- sum(scores * rel_freq)
  sigma_sq <- sum(scores^2 * rel_freq) - mu^2
  if (sigma_sq <= 0) return(sum(rel_freq[scores <= x]))
  ax <- sqrt(sigma_sq / (sigma_sq + hx^2))
  z_scores <- (x - ax * scores - (1 - ax) * mu) / (ax * hx)
  cdf_val <- sum(rel_freq * pnorm(z_scores))
  return(cdf_val)
}

#' Inverse of the Kernel Smoothed CDF
#'
#' @description
#' Finds the score `x` that corresponds to a given cumulative probability `p`.
#'
#' @param p A numeric value (probability) between 0 and 1.
#' @param scores A numeric vector of the discrete score points.
#' @param rel_freq A numeric vector of the relative frequencies for each score.
#' @param hx A numeric value for the bandwidth parameter.
#' @return The score `x` corresponding to the cumulative probability `p`.
kernel_inverse_cdf <- function(p, scores, rel_freq, hx, tol = 1e-8) {
  f_to_solve <- function(x) kernel_continu_cdf(x, scores, rel_freq, hx) - p
  search_interval <- range(scores) + c(-5 * hx, 5 * hx)
  root_result <- try(uniroot(f_to_solve, interval = search_interval, tol = tol), silent = TRUE)
  if (inherits(root_result, "try-error")) NA else root_result$root
}

# --- Optimal Bandwidth (h) Functions ---

#' First Derivative of the Kernel PDF
#' @description R translation of `KernelPdfDerivative`.
kernel_pdf_derivative <- function(x, scores, rel_freq, hx) {
  mu <- sum(scores * rel_freq)
  sigma_sq <- sum(scores^2 * rel_freq) - mu^2
  if (sigma_sq <= 0) return(0)
  ax <- sqrt(sigma_sq / (sigma_sq + hx^2))
  z_scores <- (x - ax * scores - (1 - ax) * mu) / (ax * hx)

  # Derivative of the standard normal PDF is -z * phi(z)
  derivative_val <- sum(rel_freq * (-z_scores * dnorm(z_scores)) / (ax^2 * hx^2))
  return(derivative_val)
}

#' Penalty Function for Kernel Equating Bandwidth Selection
#' @description Combines two penalty functions to find the optimal bandwidth `h`.
#' This function translates the logic of `Pen`, `Pen1`, and `Pen2`.
#' @param hx The bandwidth parameter to evaluate.
#' @param scores A numeric vector of the discrete score points.
#' @param rel_freq A numeric vector of the relative frequencies for each score.
#' @return A numeric penalty value. A lower value indicates a better fit.
penalty_function_h <- function(hx, scores, rel_freq) {
  # Penalty 1: Sum of squared differences between discrete and continuous PDF
  pdf_at_scores <- sapply(scores, kernel_continu_pdf, scores = scores, rel_freq = rel_freq, hx = hx)
  pen1 <- sum((rel_freq - pdf_at_scores)^2)

  # Penalty 2: Check for "U-shapes" in the density
  pen2 <- 0
  delta <- 1e-5 # A small offset
  for (i in seq_along(scores)) {
    deriv_left <- kernel_pdf_derivative(scores[i] - delta, scores, rel_freq, hx)
    deriv_right <- kernel_pdf_derivative(scores[i] + delta, scores, rel_freq, hx)
    # If derivative is negative on the left and positive on the right, it's a minimum (U-shape)
    if (deriv_left < 0 && deriv_right > 0) {
      pen2 <- pen2 + 1
    }
  }

  # Return combined penalty
  return(pen1 + pen2)
}

#' Find Optimal Bandwidth (h) for Kernel Equating
#' @description Finds the optimal bandwidth `h` by minimizing the penalty function.
#' This is an R translation of `Optimalh`.
#' @param scores A numeric vector of the discrete score points.
#' @param rel_freq A numeric vector of the relative frequencies for each score.
#' @param search_interval A vector of length 2 defining the search interval for `h`.
#' @return The optimal bandwidth value `h`.
find_optimal_h <- function(scores, rel_freq, search_interval = c(0.1, 2.0)) {
  optimize(
    f = penalty_function_h,
    interval = search_interval,
    scores = scores,
    rel_freq = rel_freq
  )$minimum
}

# --- Main Equating and Standard Error Functions ---

#' Perform Kernel Equipercentile Equating
#' @description This is the core kernel equating function. It takes two distributions
#' and finds the equated scores. R translation of `KernelEquate`.
#' @param scores_x, scores_y Numeric vectors of score points for Form X and Y.
#' @param rel_freq_x, rel_freq_y Numeric vectors of relative frequencies for Form X and Y.
#' @param hx, hy Optimal bandwidth parameters for Form X and Y.
#' @return A numeric vector of equated Y scores for each score on Form X.
kernel_equate <- function(scores_x, rel_freq_x, hx, scores_y, rel_freq_y, hy) {
  # 1. For each score on X, find its cumulative probability (percentile rank)
  cdf_x <- sapply(scores_x, kernel_continu_cdf, scores = scores_x, rel_freq = rel_freq_x, hx = hx)

  # 2. For each probability, find the corresponding score on Y
  equated_scores <- sapply(cdf_x, kernel_inverse_cdf, scores = scores_y, rel_freq = rel_freq_y, hx = hy)

  return(equated_scores)
}

#' Calculate the C Matrix for Delta Method Standard Errors
#' @description Computes the C matrix such that C * t(C) is the variance-covariance
#' matrix of the log-linear probabilities. R translation of `ComputeCmatrixGen`.
#' @param n_persons Total number of examinees.
#' @param design_matrix The design matrix (B) from the log-linear model.
#' @param probabilities The vector of smoothed probabilities (p) from the log-linear model.
#' @return The C matrix.
compute_c_matrix <- function(n_persons, design_matrix, probabilities) {
  B <- design_matrix
  p <- probabilities

  # Sigma = D_p - p * p'
  Sigma <- diag(p) - (p %*% t(p))

  # Var(beta) = (B' * Sigma * B)^-1 / N
  # We need a robust inverse for this
  var_beta <- MASS::ginv(t(B) %*% Sigma %*% B) / n_persons

  # C_beta * t(C_beta) = Var(beta)
  c_beta <- t(chol(var_beta))

  # C = (D_p * B) * C_beta
  C <- (diag(p) %*% B) %*% c_beta

  return(C)
}

#' Calculate Partial Derivatives for Delta Method
#' @description Calculates the partial derivative of the kernel CDF F with respect
#' to the relative frequencies r. R translation of `PartialFPartialr`.
#' @return A matrix of partial derivatives.
partial_f_partial_r <- function(scores, rel_freq, hx) {
  mu <- sum(scores * rel_freq)
  sigma_sq <- sum(scores^2 * rel_freq) - mu^2
  ax <- sqrt(sigma_sq / (sigma_sq + hx^2))

  sapply(scores, function(x_i) {
    z_scores <- (x_i - ax * scores - (1 - ax) * mu) / (ax * hx)
    pnorm(z_scores)
  })
}
