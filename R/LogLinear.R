#' Create a Design Matrix for Log-Linear Models
#'
#' @description
#' Creates a design matrix based on powers of scores for univariate or
#' bivariate log-linear models.
#'
#' @param nsu, nsv Number of score categories for variables u and v.
#' @param minu, minv, incu, incv Score range and increment for u and v.
#' @param cu, cv, cuv Degrees for polynomials and cross-products.
#' @param cpm A matrix specifying the cross-product terms.
#' @param scale Logical, if TRUE, scale the design matrix columns.
#'
#' @return A list containing the raw design matrix (`B_raw`) and the
#'   potentially scaled design matrix (`B`).
design_matrix <- function(nsu, minu, incu, nsv = 0, minv = 0, incv = 0, cu, cv = 0, cuv = 0, cpm = NULL, scale = FALSE) {

  # --- Input Validation ---
  if (nsv > 0 && incv == 0) {
    stop("Error in design_matrix: If nsv > 0, incv cannot be zero.")
  }

  ncols <- cu + cv + cuv
  ncells <- if (nsv > 0) nsu * nsv else nsu

  B_raw <- matrix(0, nrow = ncells, ncol = ncols)

  if (nsv == 0) { # Univariate case
    scores_u <- score(0:(nsu - 1), minu, incu)
    if (cu > 0) {
      for (k in 1:cu) {
        B_raw[, k] <- scores_u^k
      }
    }
  } else { # Bivariate case
    scores_u <- score(0:(nsu - 1), minu, incu)
    scores_v <- score(0:(nsv - 1), minv, incv)
    grid <- expand.grid(u = scores_u, v = scores_v)

    # u polynomials
    if (cu > 0) {
      for (k in 1:cu) B_raw[, k] <- grid$u^k
    }
    # v polynomials
    if (cv > 0) {
      for (k in 1:cv) B_raw[, cu + k] <- grid$v^k
    }
    # cross-product polynomials
    if (cuv > 0) {
      for (k in 1:cuv) {
        B_raw[, cu + cv + k] <- grid$u^cpm[k, 1] * grid$v^cpm[k, 2]
      }
    }
  }

  B <- B_raw
  if (scale) {
    B <- scale(B, center = TRUE, scale = TRUE)
  }

  return(list(B_raw = B_raw, B = B))
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

  # --- Input Validation ---
  if (inc == 0) {
    stop("Error in smooth_ull: Score increment 'inc' cannot be zero.")
  }

  design <- design_matrix(nsu = ns, minu = min, incu = inc, cu = c, scale = scale)

  iter_results <- tryCatch({
    iteration(B = design$B, B_raw = design$B_raw, nct = fd, N = n, crit = crit, max_nit = max_nit)
  }, error = function(e) {
    if (grepl("singular", e$message, ignore.case = TRUE)) {
      cli::cli_div(theme = list(span.emph = list(color = "red")))
      cli::cli_alert_danger("Log-linear model failed to converge due to a singular matrix.")
      cli::cli_alert_info("This usually means the model is over-specified for the data.")
      cli::cli_bullets(c(
        "*" = "The number of moments to fit (c) may be too high for the number of score categories (ns).",
        "*" = "Try lowering the polynomial degree.",
        "i" = "Current value: {.emph c={c}}"
      ))
      stop(e)
    } else {
      # Re-throw any other type of error
      stop(e)
    }
  })

  if (is.null(iter_results)) return(NULL) # Propagate failure

  density <- iter_results$mct / n
  crfd <- cumsum(density)
  # Assuming integer scores for perc_rank calculation
  prd <- perc_rank(x = min:(min+ns-1), min = min, max = min+ns-1, inc = inc, crfd = crfd)

  return(c(iter_results, list(density = density, crfd = crfd, prd = prd)))
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
smooth_bll <- function(n, nsu, minu, incu, nsv, minv, incv, nct, cu, cv, cuv, cpm, scale = FALSE, crit = 1e-5, max_nit = 100) {

  design <- design_matrix(nsu, minu, incu, nsv, minv, incv, cu, cv, cuv, cpm, scale)

  iter_results <- tryCatch({
    iteration(B = design$B, B_raw = design$B_raw, nct = nct, N = n, crit = crit, max_nit = max_nit)
  }, error = function(e) {
    if (grepl("singular", e$message, ignore.case = TRUE)) {
      cli::cli_div(theme = list(span.emph = list(color = "red")))
      cli::cli_alert_danger("Log-linear model failed to converge due to a singular matrix.")
      cli::cli_alert_info("This usually means the model is over-specified for the data.")
      cli::cli_bullets(c(
        "*" = "The number of moments to fit may be too high for the number of unique scores.",
        "*" = "Try lowering the polynomial degrees.",
        "i" = "Current values: {.emph cu={cu}}, {.emph cv={cv}}, {.emph cuv={cuv}}"
      ))
      stop(e)
    } else {
      # Re-throw any other type of error
      stop(e)
    }
  })

  if (is.null(iter_results)) return(NULL) # Propagate failure

  # Add the input parameters to the results list to create the full 'bivar' object
  bivar_list <- c(
    iter_results,
    list(
      n = n, nsx = nsu, minx = minu, incx = incu,
      nsv = nsv, minv = minv, incv = incv,
      cu = cu, cv = cv, cuv = cuv, cpm = cpm
    )
  )

  # --- 5. Examine the Output ---
  if (!is.null(bivar_list)) {
    cat("\n--- Smoothing Results ---\n")

    # Check 1: Beta Coefficients
    # These should be numeric values. If they are extremely large, small, or NaN,
    # it indicates a numerical issue.
    cat("\nBeta Coefficients:\n")
    print(bivar_list$Beta)

    # Check 2: Smoothed Frequencies (mct)
    # These should be non-negative. The sum should equal the original sample size (n).
    cat("\nSmoothed Frequencies (mct):\n")
    print(head(bivar_list$mct))
    cat("\nAbsolute differences between observed and predicted moments:\n")
    print(abs(round(bivar_list$n_mts - bivar_list$m_mts,6)))

    cat(sprintf("\nSum of original frequencies: %f", sum(nct_b)))
    cat(sprintf("\nSum of smoothed frequencies: %f\n", sum(bivar_list$mct)))

    # Check 3: Convergence
    # The number of iterations should be less than the maximum (40).
    cat(sprintf("\nConverged in %d iterations.\n", bivar_list$nit))

    cat(sprintf("\nLikelihood Ratio Chi-Square: %f\n", bivar_list$lrchisq))


  } else {
    cat("\n--- Function call failed --- \n")
  }

  return(bivar_list)
}
