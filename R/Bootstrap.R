# Title: R Functions for Bootstrap Standard Error Estimation
#
# Description:
# This script contains R translations of the bootstrap functions from the
# "Equating Recipes" C library, found in 'Bootstrap.c'. These functions
# are for generating bootstrap samples and calculating standard errors.

#' Generate a Non-Parametric Bootstrap Sample for a Univariate Distribution
#'
#' @description
#' Creates a bootstrap sample from an observed univariate frequency distribution.
#' This is an R replacement for the C function `Boot_USTATS`.
#'
#' @details
#' This function uses sampling from a multinomial distribution, which is an
#' efficient way to perform non-parametric bootstrapping from frequency data.
#'
#' @param n The total sample size.
#' @param freq A numeric vector of the original frequency distribution.
#'
#' @return A numeric vector representing the new bootstrap frequency distribution.
bootstrap_univariate <- function(n, freq) {
  # Calculate relative frequencies to use as probabilities for resampling
  rel_freq <- freq / sum(freq)

  # Draw one sample of size n from a multinomial distribution
  as.vector(rmultinom(n = 1, size = n, prob = rel_freq))
}

#' Generate a Non-Parametric Bootstrap Sample for a Bivariate Distribution
#'
#' @description
#' Creates a bootstrap sample from an observed bivariate frequency distribution.
#' This is an R replacement for the C function `Boot_BSTATS`.
#'
#' @param n The total sample size.
#' @param bfd A matrix representing the original bivariate frequency distribution.
#'
#' @return A matrix representing the new bootstrap bivariate frequency distribution.
bootstrap_bivariate <- function(n, bfd) {
  # Flatten the matrix into a vector of probabilities
  rel_freq <- as.vector(bfd) / sum(bfd)

  # Draw a single multinomial sample
  boot_freq_vec <- as.vector(rmultinom(n = 1, size = n, prob = rel_freq))

  # Reshape the vector back into a matrix with the original dimensions
  matrix(boot_freq_vec, nrow = nrow(bfd), ncol = ncol(bfd))
}

#' Generate a Parametric Bootstrap Sample (Univariate)
#'
#' @description
#' Creates a parametric bootstrap sample from a smoothed univariate distribution.
#' This is an R replacement for `Parametric_boot_univ_BB` and
#' `Parametric_boot_univ_ULL`.
#'
#' @param n The total sample size for the bootstrap sample.
#' @param smoothed_dist A numeric vector of smoothed probabilities (must sum to 1).
#'
#' @return A numeric vector representing the new bootstrap frequency distribution.
bootstrap_parametric_univ <- function(n, smoothed_dist) {
  # Probabilities are taken directly from the smoothed distribution
  as.vector(rmultinom(n = 1, size = n, prob = smoothed_dist))
}

#' Generate a Parametric Bootstrap Sample (Bivariate)
#'
#' @description
#' Creates a parametric bootstrap sample from a smoothed bivariate distribution.
#' This is an R replacement for `Parametric_boot_biv`.
#'
#' @param n The total sample size.
#' @param smoothed_bfd A matrix of smoothed bivariate probabilities.
#'
#' @return A matrix representing the new bootstrap bivariate frequency distribution.
bootstrap_parametric_biv <- function(n, smoothed_bfd) {
  # Flatten the smoothed matrix into a vector of probabilities
  rel_freq <- as.vector(smoothed_bfd)

  # Draw a single multinomial sample
  boot_freq_vec <- as.vector(rmultinom(n = 1, size = n, prob = rel_freq))

  # Reshape the vector back into a matrix
  matrix(boot_freq_vec, nrow = nrow(smoothed_bfd), ncol = ncol(smoothed_bfd))
}


#' Calculate Bootstrap Statistics from Replications
#'
#' @description
#' Calculates the mean, standard deviation (bootstrap SE), percentile-based
#' confidence intervals, and overall SE from a matrix of bootstrap results.
#' This function replaces the C functions `Boot_accumulate_eraw`, `Boot_se_eraw`,
#' `Boot_accumulate_ess`, and `Boot_se_ess`.
#'
#' @details
#' In R, it is more idiomatic to collect all results into a matrix first and
#' then compute summary statistics, rather than accumulating sums in a loop.
#'
#' @param results_matrix A matrix where each column is the result from one
#'   bootstrap replication and each row corresponds to a score point.
#' @param original_freq (Optional) A vector of original frequencies, used to
#'   calculate the overall weighted standard error.
#'
#' @return A list containing `mean_scores`, `se_scores`, `overall_se`, and
#'   the 95% confidence interval bounds `ci_95_low` and `ci_95_high`.
calculate_bootstrap_se <- function(results_matrix, original_freq = NULL) {
  # Calculate the mean and standard deviation across replications for each score point
  mean_scores <- apply(results_matrix, 1, mean)
  ci_95_low <- apply(results_matrix, 1, quantile, .025)
  ci_95_high <- apply(results_matrix, 1, quantile, .975)
  se_scores <- apply(results_matrix, 1, sd)

  overall_se <- NA
  if (!is.null(original_freq)) {
    # Calculate the overall bootstrap SE, weighted by the original frequencies
    variances <- se_scores^2
    n <- sum(original_freq)
    if (n > 0) {
      overall_se <- sqrt(sum(original_freq * variances) / n)
    }
  }

  return(list(
    mean_scores = mean_scores,
    se_scores = se_scores,
    overall_se = overall_se,
    ci_95_low = ci_95_low,
    ci_95_high = ci_95_high
  ))
}
