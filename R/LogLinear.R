#' Create a Design Matrix for Log-Linear Models
#'
#' @description
#' Creates a design matrix and returns the scaling attributes (center and scale)
#' used during its creation, which are necessary for unscaling the results.
#' @return A list containing the raw design matrix (`B_raw`), the scaled design
#'   matrix (`B`), and the scaling attributes (`B_center`, `B_scale`).
design_matrix <- function(nsu, minu, incu, nsv = 0, minv = 0, incv = 0, cu, cv = 0, cuv = 0, cpm = NULL, scale = FALSE) {
  ncols <- cu + cv + cuv
  ncells <- if (nsv > 0) nsu * nsv else nsu
  B_raw <- matrix(0, nrow = ncells, ncol = ncols)

  score <- function(loc, min, inc) min + loc * inc

  if (nsv == 0) { # Univariate
    scores_u <- score(0:(nsu - 1), minu, incu)
    if (cu > 0) {
      for (k in 1:cu) B_raw[, k] <- scores_u^k
    }
  } else { # Bivariate
    scores_u <- score(0:(nsu - 1), minu, incu)
    scores_v <- score(0:(nsv - 1), minv, incv)
    grid <- expand.grid(u = scores_u, v = scores_v)
    if (cu > 0) for (k in 1:cu) B_raw[, k] <- grid$u^k
    if (cv > 0) for (k in 1:cv) B_raw[, cu + k] <- grid$v^k
    if (cuv > 0) for (k in 1:cuv) B_raw[, cu + cv + k] <- grid$u^cpm[k, 1] * grid$v^cpm[k, 2]
  }

  B <- B_raw
  B_center <- rep(0, ncols)
  B_scale <- rep(1, ncols)

  if (scale) {
    B_scaled <- scale(B_raw, center = TRUE, scale = TRUE)
    B <- as.matrix(B_scaled)
    B_center <- attr(B_scaled, "scaled:center")
    B_scale <- attr(B_scaled, "scaled:scale")
    B_scale[B_scale == 0] <- 1 # Avoid division by zero
  }
  return(list(B_raw = B_raw, B = B, B_center = B_center, B_scale = B_scale))
}

#' Convert Scaled Beta Coefficients to Unscaled Metric
#'
#' @description
#' Transforms the log-linear model's Beta coefficients and intercept (ap)
#' from a scaled metric back to the original, unscaled metric.
#' @return A list containing the unscaled Beta coefficients and intercept.
unscale_beta <- function(Beta_scaled, ap_scaled, B_center, B_scale) {
  Beta_unscaled <- Beta_scaled / B_scale
  ap_unscaled <- ap_scaled - sum(Beta_scaled * B_center / B_scale)
  return(list(Beta_unscaled = Beta_unscaled, ap_unscaled = ap_unscaled))
}

#' Get First Derivative of Log-Likelihood
#'
#' @description Computes B' * (n - m).
get_Btnm <- function(B, n, m) {
  return(t(B) %*% (n - m))
}

#' Get Negative Second Derivative of Log-Likelihood
#'
#' @description Computes B' * S_m * B.
get_BtSmB <- function(B, m, N) {
  Bm <- t(B) %*% m
  return(t(B) %*% diag(m) %*% B - (Bm %*% t(Bm)) / N)
}

#' Get Initial Beta Coefficients
#'
#' @description Computes initial estimates for the log-linear model coefficients.
get_Beta0 <- function(B, n, N) {
  ns <- length(n)
  a <- 0.8 * n + 0.2 * N / ns

  BtSaB <- get_BtSmB(B, a, N)

  aloga <- sum(a * log(a))
  Baloga <- t(B) %*% (a * log(a))
  Ba <- t(B) %*% a
  BtSaloga <- Baloga - Ba * aloga / N

  # Solve the system BtSaB * Beta = BtSaloga
  beta0 <- er_lubksb(a = BtSaB, b = BtSaloga)
  return(as.vector(beta0))
}

#' Get Fitted Frequencies (mct)
#'
#' @description Computes the fitted frequencies from the Beta coefficients.
get_mct <- function(B, Beta, N) {
  BBeta <- B %*% Beta
  ap <- log(N) - log(sum(exp(BBeta)))
  m <- exp(ap + BBeta)
  return(list(mct = as.vector(m), ap = ap))
}

#' Check Convergence Criterion
#'
#' @description Checks if the moments of the observed and fitted distributions have converged.
crit_mts <- function(n_mts, m_mts, crit) {
  all(abs(n_mts - m_mts) <= crit)
}

#' Get Log-Linear Moments
#'
#' @description Calculates moments based on the design matrix and frequencies.
get_LLmoments <- function(B, B_raw, f, N) {
  rel_freq <- f / N
  mts <- t(B) %*% rel_freq

  # Simplified version for raw moments for this example
  mts_raw <- t(B_raw) %*% rel_freq

  return(list(mts = as.vector(mts), mts_raw = as.vector(mts_raw)))
}

#' Main Iterative Algorithm for Log-Linear Fitting
#'
#' @description The Newton-Raphson algorithm to fit the log-linear model.
iteration <- function(B, B_raw, nct, N, max_nit = 100, crit = 1e-5) {

  nc <- ncol(B)

  # Get initial beta
  Beta <- get_Beta0(B, nct, N)
  converged <- FALSE # Initialize convergence flag


  for (nit in 1:max_nit) {
    # Get fitted frequencies
    m_info <- get_mct(B, Beta, N)
    mct <- m_info$mct

    # Get moments
    n_moments <- get_LLmoments(B, B_raw, nct, N)
    m_moments <- get_LLmoments(B, B_raw, mct, N)

    # Check for convergence
    if (crit_mts(n_moments$mts, m_moments$mts, crit)) {
      converged <- TRUE
      break
    }

    # Update Beta if not converged
    BtSmB <- get_BtSmB(B, mct, N)
    Btnm <- get_Btnm(B, nct, mct)

    delta <- er_lubksb(BtSmB, Btnm)
    Beta <- Beta + delta
  }

  if (nit == max_nit) warning("Max iterations reached without convergence.")

  lrchisq <- 2 * sum(nct[nct > 0] * log(nct[nct > 0] / mct[nct > 0]))

  return(list(
    Beta = Beta, mct = mct, nit = nit, lrchisq = lrchisq,
    n_mts = n_moments$mts, m_mts = m_moments$mts, ap = m_info$ap, converged = converged
  ))
}

#' Perform Univariate Log-Linear Smoothing
#'
#' @description
#' This is a high-level wrapper function that performs univariate log-linear
#' smoothing on a frequency distribution. It creates the design matrix, runs
#' the iterative fitting algorithm, and returns the smoothed distribution and
#' related statistics. This is an R translation of `Smooth_ULL`.
#'
#' @param n An integer, the total number of examinees.
#' @param ns An integer, the number of score categories.
#' @param min A numeric value, the minimum score on the scale.
#' @param inc A numeric value, the increment between scores.
#' @param fd A numeric vector of the observed score frequencies.
#' @param c An integer, the degree of the polynomial to fit (e.g., c=3 for a
#'   cubic model). This determines the number of moments to be matched.
#' @param scale A logical value. If TRUE, the columns of the design matrix
#'   are centered and scaled, which can improve numerical stability.
#' @param crit A numeric value, the convergence criterion for the moments.
#'
#' @return A list containing the full results of the smoothing process,
#'   including the Beta coefficients, fitted frequencies (`mct`), number of
#'   iterations, chi-square statistics, moments, and the final smoothed
#'   density, CDF, and percentile ranks.
#' @seealso \code{\link{iteration}}, \code{\link{design_matrix}}
smooth_ull <- function(n, ns, min, inc, fd, c, scale = FALSE, crit = 1e-5, max_nit = 100) {
  # Always use scaling internally for numerical stability
  design <- design_matrix(nsu = ns, minu = min, incu = inc, cu = c, scale = TRUE)

  iter_results <- iteration(B = design$B, B_raw = design$B_raw, nct = fd, N = n, crit = crit, max_nit = max_nit)

  unscaled_params <- unscale_beta(
    Beta_scaled = iter_results$Beta,
    ap_scaled = iter_results$ap,
    B_center = design$B_center,
    B_scale = design$B_scale
  )

  iter_results$Beta <- unscaled_params$Beta_unscaled
  iter_results$ap <- unscaled_params$ap_unscaled

  density <- iter_results$mct / n
  crfd <- cumsum(density)

  return(c(iter_results, list(density = density, crfd = crfd)))
}


#' Perform Bivariate Log-Linear Smoothing
#'
#' @description
#' This is a high-level wrapper function that performs bivariate log-linear
#' smoothing on a two-way frequency distribution. It creates the design matrix,
#' runs the iterative fitting algorithm, and returns the results. This is an R
#' translation of `Smooth_BLL`.
#'
#' @param n An integer, the total number of examinees.
#' @param nsu An integer, the number of score categories for the first variable (u).
#' @param minu A numeric value, the minimum score for variable u.
#' @param incu A numeric value, the increment for variable u.
#' @param nsv An integer, the number of score categories for the second variable (v).
#' @param minv A numeric value, the minimum score for variable v.
#' @param incv A numeric value, the increment for variable v.
#' @param nct A numeric vector representing the flattened bivariate frequency
#'   table. The vector should be created by reading the table row-by-row
#'   (i.e., `as.vector(t(bivariate_matrix))`).
#' @param cu An integer, the degree of the polynomial for variable u.
#' @param cv An integer, the degree of the polynomial for variable v.
#' @param cuv An integer, the number of cross-product terms in the model.
#' @param cpm A matrix with `cuv` rows and 2 columns, where each row specifies
#'   the powers `(i, j)` for a cross-product term `u^i * v^j`.
#' @param scale A logical value. If TRUE, the columns of the design matrix
#'   are centered and scaled.
#' @param crit A numeric value, the convergence criterion for the moments.
#'
#' @return A list containing the full results of the smoothing process,
#'   including the Beta coefficients, fitted frequencies (`mct`), number of
#'   iterations, chi-square statistics, and moments. This list is the
#'   `bivar` object required by the `equate_cll` function.
#' @seealso \code{\link{iteration}}, \code{\link{design_matrix}}, \code{\link{equate_cll}}

#' Perform Bivariate Log-Linear Smoothing
smooth_bll <- function(n, nsu, minu, incu, nsv, minv, incv, nct, cu, cv, cuv, cpm, scale = FALSE, crit = 1e-5, max_nit = 100) {
  # Step 1: Always create a scaled design matrix for internal use
  design <- design_matrix(nsu, minu, incu, nsv, minv, incv, cu, cv, cuv, cpm, scale = TRUE)

  # Step 2: Run the iteration with the scaled design matrix
  iter_results <- tryCatch({
    iteration(B = design$B, B_raw = design$B_raw, nct = nct, N = n, crit = crit, max_nit = max_nit)
  }, error = function(e) {
    if (grepl("singular", e$message, ignore.case = TRUE)) {
      cli::cli_div(theme = list(span.emph = list(color = "red")))
      cli::cli_alert_danger("Log-linear model failed to converge due to a singular matrix.")
      cli::cli_alert_info("This can happen if the model is over-specified.")
      cli::cli_bullets(c("*" = "Try lowering the polynomial degrees."))
      stop(e)
    } else {
      stop(e)
    }
  })

  if (is.null(iter_results)) return(NULL)

  # Step 3: Conditionally un-scale the coefficients based on the `scale` argument
  if (scale == FALSE) {
    # If the user wants unscaled coefficients, transform them back
    unscaled_params <- unscale_beta(
      Beta_scaled = iter_results$Beta,
      ap_scaled = iter_results$ap,
      B_center = design$B_center,
      B_scale = design$B_scale
    )
    # Overwrite the scaled results with the unscaled ones
    iter_results$Beta <- unscaled_params$Beta_unscaled
    iter_results$ap <- unscaled_params$ap_unscaled
  }
  # If scale == TRUE, do nothing. iter_results already contains the scaled coefficients.

  # Step 4: Assemble the final `bivar` object
  bivar_list <- c(
    iter_results,
    list(
      n = n, nsx = nsu, minx = minu, incx = incu,
      nsv = nsv, minv = minv, incv = incv,
      cu = cu, cv = cv, cuv = cuv, cpm = cpm,
      # Add a flag to indicate the state of the returned model
      scaled_model = scale
    )
  )

  # --- Diagnostic Printing (Merged from your version) ---
  if (!is.null(bivar_list)) {
    cat("\n--- Smoothing Results ---\n")
    cat("\nBeta Coefficients:\n")
    print(bivar_list$Beta)
    cat("\nSmoothed Frequencies (mct):\n")
    print(head(bivar_list$mct))
    cat("\nAbsolute differences between observed and predicted moments:\n")
    print(abs(round(bivar_list$n_mts - bivar_list$m_mts, 6)))
    cat(sprintf("\nSum of original frequencies: %f", sum(nct)))
    cat(sprintf("\nSum of smoothed frequencies: %f\n", sum(bivar_list$mct)))
    cat(sprintf("\nConverged in %d iterations.\n", bivar_list$nit))
    cat(sprintf("\nLikelihood Ratio Chi-Square: %f\n", bivar_list$lrchisq))
  } else {
    cat("\n--- Function call failed --- \n")
  }

  return(bivar_list)
}
