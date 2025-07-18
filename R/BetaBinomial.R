# Title: R Functions for Beta-Binomial Smoothing
#
# Description:
# This script contains R translations of the functions from 'BetaBinomial.c'.
# These functions implement the beta-binomial smoothing method for score
# distributions, as described by Hanson (1991) and referenced in
# Kolen & Brennan (2004).

# Source the file containing previously converted utility functions.
# Note: You may need to change this file path.
# source("equating_recipes_utils.R")

#' Calculate Lord's k from KR20 Reliability
#'
#' @description
#' Calculates the value of Lord's k for the compound binomial error model,
#' given a test's reliability (typically KR-20), number of items, and raw
#' score moments. This is an R translation of `CalcLordk`.
#'
#' @param kr20 A numeric value for the reliability coefficient (e.g., KR-20).
#' @param n_items An integer for the number of items on the test.
#' @param rmoment A numeric vector of the first two raw score moments
#'   (mean and standard deviation).
#'
#' @return The numeric value of Lord's k.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
calc_lord_k <- function(kr20, n_items, rmoment) {
  dn <- as.double(n_items)
  mnm <- rmoment[1] * (dn - rmoment[1])
  varr <- rmoment[2]^2

  varp <- mnm / (dn^2)
  varp <- varp - (varr / dn) * (1.0 - ((dn - 1.0) / dn) * kr20)

  k <- mnm - varr - dn * varp
  if (k == 0) return(0) # Avoid division by zero
  k <- (dn * dn * (dn - 1.0) * varp) / (2.0 * k)
  return(k)
}

#' Compute True Score Moments from Raw Score Moments
#'
#' @description
#' Derives the first four moments of the unobserved true-score distribution
#' from the moments of the observed raw-score distribution, based on Lord's
#' compound binomial error model. This is an R translation of `BetaMoments`.
#'
#' @param n An integer, the number of items.
#' @param k A numeric value for Lord's k parameter. Use k=0 for the simple
#'   binomial error model.
#' @param rmoment A numeric vector of the first four raw score moments
#'   (mean, sd, skewness, kurtosis).
#'
#' @return A list containing two numeric vectors:
#'   \item{tmoment}{The first four central moments (mean, sd, skew, kurtosis)
#'     of the estimated true-score distribution.}
#'   \item{nctmoment}{The first four non-central moments of the estimated
#'     true-score distribution.}
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
beta_moments <- function(n, k, rmoment) {
  # Convert raw score central moments to non-central moments
  mu1 <- rmoment[1]
  mu2 <- rmoment[2]^2 + mu1^2
  mu3 <- rmoment[3] * rmoment[2]^3 + 3*mu1*rmoment[2]^2 + mu1^3
  mu4 <- rmoment[4] * rmoment[2]^4 + 4*rmoment[3]*rmoment[2]^3*mu1 + 6*mu1^2*rmoment[2]^2 + mu1^4

  # Convert non-central moments to factorial moments
  f1 <- mu1
  f2 <- mu2 - mu1
  f3 <- mu3 - 3*mu2 + 2*mu1
  f4 <- mu4 - 6*mu3 + 11*mu2 - 6*mu1

  # Estimate true-score non-central moments from factorial moments
  n_dbl <- as.double(n)
  nctmoment <- numeric(4)
  nctmoment[1] <- f1 / n_dbl
  nctmoment[2] <- (f2 + 2*k*nctmoment[1]) / (n_dbl * (n_dbl - 1) + 2*k)
  nctmoment[3] <- (f3 + 6*k*nctmoment[2]) / (n_dbl * (n_dbl - 1) * (n_dbl - 2) + 6*k)
  nctmoment[4] <- (f4 + 12*k*nctmoment[3]) / (n_dbl * (n_dbl - 1) * (n_dbl - 2) * (n_dbl-3) + 12*k)

  # Convert true-score non-central moments back to central moments
  tmoment <- numeric(4)
  tmoment[1] <- nctmoment[1] * n_dbl
  t_var <- nctmoment[2] * n_dbl^2 - tmoment[1]^2

  if (t_var > 0) {
    tmoment[2] <- sqrt(t_var)
    t_cm3 <- nctmoment[3]*n_dbl^3 - 3*nctmoment[2]*n_dbl^2*tmoment[1] + 2*tmoment[1]^3
    tmoment[3] <- t_cm3 / (t_var^(3/2))
    t_cm4 <- nctmoment[4]*n_dbl^4 - 4*nctmoment[3]*n_dbl^3*tmoment[1] + 6*nctmoment[2]*n_dbl^2*tmoment[1]^2 - 3*tmoment[1]^4
    tmoment[4] <- t_cm4 / (t_var^2)
  } else {
    tmoment[2:4] <- 0
  }

  return(list(tmoment = tmoment, nctmoment = nctmoment))
}


#' Calculate Kurtosis from Beta Parameters
#'
#' @description
#' Calculates the kurtosis of a beta distribution given its two shape parameters.
#' This is an R translation of `CalcKurt`.
#'
#' @param alpha The first shape parameter (alpha > 0).
#' @param beta The second shape parameter (beta > 0).
#'
#' @return The numeric value of the kurtosis.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
calc_kurt <- function(alpha, beta) {
  a <- alpha; b <- beta
  k1 <- 3.0 * (a + b + 1.0)
  k2 <- 2.0 * (a + b)^2
  k3 <- a * b
  k4 <- a + b - 6.0
  k5 <- a + b + 2.0
  k6 <- a + b + 3.0
  kurt <- (k1 * (k2 + k3 * k4)) / (k3 * k5 * k6)
  return(kurt)
}

#' Estimate Negative Hypergeometric Parameters (2-Parameter Beta)
#'
#' @description
#' Estimates the alpha and beta shape parameters for a 2-parameter beta
#' distribution (or a scaled 4-parameter beta) from true score moments.
#' This is an R translation of `EstNegHypGeo`.
#'
#' @param n_items The number of items.
#' @param tmoment A numeric vector of the first two true score moments (mean, sd).
#' @param lower The lower bound of the beta distribution (e.g., 0).
#' @param upper The upper bound of the beta distribution (e.g., 1).
#'
#' @return A named vector with `alpha` and `beta` parameters, or `NULL` if
#'   estimation fails.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
est_neg_hyp_geo <- function(n_items, tmoment, lower, upper) {
  maxmin <- n_items * (upper - lower)
  if(maxmin == 0) return(NULL)
  smean <- (tmoment[1] - n_items * lower) / maxmin
  svar <- (tmoment[2]^2) / (maxmin^2)
  if(svar <= 0) return(NULL)
  alpha <- smean^2 * (1.0 - smean) / svar - smean
  beta <- alpha * (1.0 - smean) / smean
  if (is.finite(alpha) && is.finite(beta) && alpha > 0 && beta > 0) {
    return(c(alpha = alpha, beta = beta))
  } else {
    return(NULL)
  }
}

#' Estimate 4-parameter Beta via Method of Moments
#'
#' @description
#' Estimates the four parameters of a beta distribution by matching the first
#' four moments of the true score distribution. This is an R translation of
#' `CalcBetaParaMM`.
#'
#' @param n_items The number of items.
#' @param tmoment A numeric vector of the first four true score moments.
#'
#' @return A named vector with `alpha`, `beta`, `lower`, and `upper` parameters,
#'   or `NULL` if estimation fails.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
calc_beta_params_mm <- function(n_items, tmoment) {
  b2 <- tmoment[4]; b1 <- tmoment[3]^2
  denom <- 6.0 + 3.0 * b1 - 2.0 * b2
  if (denom == 0) return(NULL)
  r <- 6.0 * (b2 - b1 - 1.0) / denom
  rad_denom <- (r + 2.0) * (r + 3.0) * b2 - 3.0 * (r + 1.0) * (r - 6.0)
  if (rad_denom == 0) return(NULL)
  rad <- (24.0 * (r + 1.0)) / rad_denom
  if (rad > 1.0 || rad < 0) return(NULL)
  rad <- sqrt(1.0 - rad); r <- r / 2.0
  a <- if (tmoment[3] < 0.0) r * (1.0 + rad) else r * (1.0 - rad)
  b <- if (tmoment[3] < 0.0) r * (1.0 - rad) else r * (1.0 + rad)
  if (a <= 0.0 || b <= 0.0) return(NULL)
  bma <- (a + b) * sqrt(a + b + 1) * tmoment[2] / sqrt(a * b)
  l <- -bma * (a / (a + b)) + tmoment[1]
  u <- bma + l
  lower <- l / n_items; upper <- u / n_items
  if (lower < 0.0 || upper > 1.0) return(NULL)
  return(c(alpha = a, beta = b, lower = lower, upper = upper))
}

#' Estimate 4-parameter Beta via Least Squares
#'
#' @description
#' Finds the beta distribution parameters that match the first three true score
#' moments while minimizing the squared difference between the observed and
#' fitted kurtosis. This is an R replacement for `CalcBetaParaLS`.
#'
#' @param n_items The number of items.
#' @param tmoment A vector of the first four true score moments.
#' @param nctmoment A vector of the first four non-central true score moments.
#'
#' @return A named vector with the four beta parameters, or `NULL` on failure.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
calc_beta_params_ls <- function(n_items, tmoment, nctmoment) {

  find_upper <- function(lower, nctm) {
    m1 <- nctm[1]; m2 <- nctm[2]; m3 <- nctm[3]
    num <- m1^2 * (m2 * lower - 2*m3) + m2^2 * (m1 - 2*lower) + m3 * (m2 + m1*lower)
    den <- m1^2 * (2*m1*lower - m2) + m2 * (2*m2 - 3*m1*lower) + m3 * (lower - m1)
    if (den == 0) return(NA)
    return(num / den)
  }

  kurt_diff_sq <- function(lower, n_items, tmoment, nctmoment) {
    upper <- find_upper(lower, nctmoment)
    if (is.na(upper) || upper <= lower || upper > 1.0) return(1e10)

    params <- est_neg_hyp_geo(n_items, tmoment, lower, upper)
    if (is.null(params)) return(1e10)

    kurt_fit <- calc_kurt(params["alpha"], params["beta"])
    return((tmoment[4] - kurt_fit)^2)
  }

  optim_res <- try(optimize(kurt_diff_sq, interval = c(0, 1), n_items = n_items, tmoment = tmoment, nctmoment = nctmoment), silent=TRUE)
  if(inherits(optim_res, "try-error")) return(NULL)

  best_lower <- optim_res$minimum
  best_upper <- find_upper(best_lower, nctmoment)

  params <- est_neg_hyp_geo(n_items, tmoment, best_lower, best_upper)
  if (is.null(params)) return(NULL)

  return(c(params, lower = best_lower, upper = best_upper))
}


#' Estimate Beta Distribution Parameters
#'
#' @description
#' Main dispatcher for parameter estimation. It tries to fit the specified number
#' of parameters, with fallbacks to simpler models if estimation fails. This is
#' an R translation of `CalcBetaPara`.
#'
#' @param n_items The number of items.
#' @param tmoment A vector of the first four true score moments.
#' @param nctmoment A vector of the first four non-central true score moments.
#' @param nparm The desired number of parameters to fit (2 or 4).
#'
#' @return A named vector with the final estimated parameters (`alpha`, `beta`,
#'   `lower`, `upper`) and `moments_fit`, indicating which method succeeded.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
calc_beta_params <- function(n_items, tmoment, nctmoment, nparm) {
  if (nparm == 4 && tmoment[2] > 0) {
    params <- calc_beta_params_mm(n_items, tmoment)
    if (!is.null(params)) return(c(params, moments_fit = 4))

    params <- calc_beta_params_ls(n_items, tmoment, nctmoment)
    if (!is.null(params)) return(c(params, moments_fit = 3))
  }

  params <- est_neg_hyp_geo(n_items, tmoment, 0, 1)
  if (!is.null(params)) return(c(params, lower = 0, upper = 1, moments_fit = 2))

  alpha <- 1.0
  beta <- (n_items / tmoment[1]) - alpha
  if (beta > 0) return(c(alpha = alpha, beta = beta, lower = 0, upper = 1, moments_fit = 1))

  return(c(alpha = 1, beta = 1, lower = 0, upper = 1, moments_fit = 0))
}

#' Calculate Observed Density from Beta Parameters
#'
#' @description
#' Calculates the expected beta-binomial frequency distribution using the
#' `betafunctions` package, which can handle the 4-parameter model.
#'
#' @param n_items The number of items.
#' @param n_persons The total sample size.
#' @param beta_params A named vector with beta distribution parameters:
#'   `alpha`, `beta`, `lower`, `upper`.
#'
#' @return A vector of expected frequencies for scores 0 to n_items.
obs_density <- function(n_items, n_persons, beta_params) {
  if (!requireNamespace("betafunctions", quietly = TRUE)) {
    stop("Package 'betafunctions' is needed for this function. Please install it.", call. = FALSE)
  }
  probs <- betafunctions::dBetaBinom(
    x = 0:n_items,
    N = n_items,
    l = beta_params["lower"],
    u = beta_params["upper"],
    alpha = beta_params["alpha"],
    beta = beta_params["beta"]
  )
  return(probs * n_persons)
}

#' Adjust Density for Compound Binomial Error
#'
#' @description
#' Adjusts an observed score distribution for compound binomial error using
#' Lord's k. This is an R translation of `ObsDenK`.
#'
#' @param n_items The number of items.
#' @param k The numeric value of Lord's k.
#' @param scounts A numeric vector of fitted frequencies from `obs_density`.
#'
#' @return A numeric vector of the adjusted frequencies.
#' @author B. A. Hanson (Original C code), Google's Gemini (R translation)
obs_den_k <- function(n_items, k, scounts) {
  dn <- as.double(n_items)
  kn2 <- k / (dn * (dn - 1.0))

  px <- (dn - 1.0) * scounts[2]
  scounts[1] <- scounts[1] - kn2 * px
  pxm1 <- 0.0

  if (n_items > 1) {
    for (s in 2:n_items) {
      dsp1 <- as.double(s)
      pxp1 <- dsp1 * (dn - dsp1) * scounts[s + 1]
      d2 <- pxm1 - 2.0 * px + pxp1
      scounts[s] <- scounts[s] - kn2 * d2
      pxm1 <- px
      px <- pxp1
    }
  }

  scounts[n_items + 1] <- scounts[n_items + 1] - kn2 * pxm1
  return(scounts)
}

#' Likelihood Ratio Chi-Square
#'
#' @description
#' Computes the likelihood ratio chi-squared statistic to assess goodness-of-fit.
#' This is an R translation of `LRChiSqr`.
#'
#' @param raw_counts A vector of observed frequencies.
#' @param fit_counts A vector of fitted frequencies from the model.
#'
#' @return The numeric chi-square value.
lr_chi_sqr <- function(raw_counts, fit_counts) {
  fit_counts[fit_counts == 0] <- 1e-9
  terms <- ifelse(raw_counts > 0, raw_counts * log(raw_counts / fit_counts), 0)
  return(2 * sum(terms))
}

#' Pearson Chi-Square
#'
#' @description
#' Computes the Pearson chi-squared statistic to assess goodness-of-fit.
#' This is an R translation of `PChiSqr`.
#'
#' @param raw_counts A vector of observed frequencies.
#' @param fit_counts A vector of fitted frequencies from the model.
#'
#' @return The numeric chi-square value.
p_chi_sqr <- function(raw_counts, fit_counts) {
  fit_counts[fit_counts == 0] <- 1e-9
  return(sum((raw_counts - fit_counts)^2 / fit_counts))
}

#' Perform Beta-Binomial Smoothing
#'
#' @description
#' This is the main wrapper function that orchestrates the entire beta-binomial
#' smoothing process. It estimates true score moments, fits a beta distribution
#' to those moments, and then calculates the smoothed observed score distribution.
#'
#' @param n_persons The total number of examinees in the sample.
#' @param n_items The number of items on the test.
#' @param freq A numeric vector of the observed raw score frequencies.
#' @param rmoment A numeric vector of the first four raw score moments.
#' @param nparm The number of parameters for the beta distribution (2 or 4).
#' @param rel The reliability of the test (e.g., KR-20), used to calculate
#'   Lord's k for the compound binomial model. Set to 0 for the simple model.
#'
#' @return A list containing all results from the smoothing process, including
#'   the smoothed density, beta parameters, fitted moments, and chi-square statistics.
#' @author B. A. Hanson & R. L. Brennan (Original C code), Google's Gemini (R translation)
smooth_bb <- function(n_persons, n_items, freq, rmoment, nparm, rel) {

  lord_k <- if (rel > 0) calc_lord_k(rel, n_items, rmoment) else 0
  if (is.na(lord_k) || lord_k < 0) lord_k <- 0

  moments <- beta_moments(n_items, lord_k, rmoment)
  tmoment <- moments$tmoment
  nctmoment <- moments$nctmoment

  params_fit <- calc_beta_params(n_items, tmoment, nctmoment, nparm)
  beta_params <- params_fit[1:4]

  fitted_counts <- obs_density(n_items, n_persons, beta_params)

  if (lord_k > 0) {
    fitted_counts <- obs_den_k(n_items, lord_k, fitted_counts)
  }

  fitted_counts[fitted_counts < 0] <- 0

  smoothed_density <- fitted_counts / sum(fitted_counts)

  fitted_moments <- get_moments(scores = 0:n_items, rel_freq = smoothed_density)

  lr_chisq <- lr_chi_sqr(freq, fitted_counts)
  p_chisq <- p_chi_sqr(freq, fitted_counts)

  return(list(
    density = smoothed_density,
    crfd = cumsum(smoothed_density),
    prd = perc_rank(x = 0:n_items, min = 0, max = n_items, inc = 1, crfd = cumsum(smoothed_density)),
    beta_params = beta_params,
    fitted_moments = fitted_moments,
    true_moments = tmoment,
    lr_chisq = lr_chisq,
    p_chisq = p_chisq,
    moments_fit = params_fit["moments_fit"]
  ))
}
