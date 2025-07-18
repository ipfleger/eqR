
#' Calculate Synthetic Population Moments for Linear Equating (CINEG)
#'
#' @description
#' A helper function for CINEG linear equating. It calculates the means and
#' standard deviations of the synthetic population, as well as the slope and
#' intercept of the linear equating function. This is an R translation of the
#' `CI_LinObsEq` function.
linear_equate_observed_ci <- function(mnx1, sdx1, mnv1, sdv1, mny2, sdy2, mnv2, sdv2, w1, mean_only, gamma1, gamma2) {
  w2 <- 1 - w1

  musx <- mnx1 - w2 * gamma1 * (mnv1 - mnv2)
  musy <- mny2 + w1 * gamma2 * (mnv1 - mnv2)

  if (!mean_only) {
    varsx <- sdx1^2 - w2 * gamma1^2 * (sdv1^2 - sdv2^2) + w1 * w2 * gamma1^2 * (mnv1 - mnv2)^2
    varsy <- sdy2^2 + w1 * gamma2^2 * (sdv1^2 - sdv2^2) + w1 * w2 * gamma2^2 * (mnv1 - mnv2)^2

    if (varsx <= 0) stop("Synthetic variance for X is non-positive.")

    a <- sqrt(varsy / varsx)
    ssx <- sqrt(varsx)
    ssy <- sqrt(varsy)
  } else {
    a <- 1.0
    ssx <- NA
    ssy <- NA
  }

  b <- musy - a * musx

  return(list(msx = musx, msy = musy, ssx = ssx, ssy = ssy, a = a, b = b))
}

#' Perform Linear Equating for the Common-Item Non-Equivalent Groups (CINEG) Design
#'
#' @description
#' Computes linear equating for the CINEG design using specified methods:
#' Tucker, Levine Observed Score, Levine True Score, and/or Chained Linear. This is
#' an R translation of the `CI_LinEq` function.
#'
#' @param mnx1,sdx1 Mean and SD of form X for population 1.
#' @param mnv1,sdv1 Mean and SD of anchor test V for population 1.
#' @param covxv1 Covariance of X and V for population 1.
#' @param mny2,sdy2 Mean and SD of form Y for population 2.
#' @param mnv2,sdv2 Mean and SD of anchor test V for population 2.
#' @param covyv2 Covariance of Y and V for population 2.
#' @param w1 Weight for population 1 in the synthetic population.
#' @param anchor Is the anchor internal (1) or external (0)?
#' @param mean_only A logical. If TRUE, mean equating is performed (slope fixed to 1).
#'   If FALSE (default), linear equating is performed.
#' @param type A character vector specifying which methods to run. Can include
#'   `"tucker"`, `"levine_observed"`, `"levine_true"`, `"chained"`. The default
#'   `"all"` runs all four methods.
#' @param min_x,max_x,inc_x The score scale for form X.
#'
#' @return A list containing a `summary` data frame with coefficients for each
#'   method run, and an `equated_scores` data frame with the equated scores.
linear_equate_ci <- function(mnx1, sdx1, mnv1, sdv1, covxv1, mny2, sdy2, mnv2, sdv2, covyv2, w1, anchor = FALSE, mean_only = FALSE, type = "all", min_x, max_x, inc_x) {

  # if (w1 < 0) {
  #   if (is.null(n1) || is.null(n2)) {
  #     stop("If w1 is negative, n1 and n2 must be provided to calculate proportional weights.")
  #   }
  #   w1 <- n1 / (n1 + n2)
  # }

  # Determine which methods to run
  all_methods <- c("tucker", "levine_observed", "levine_true", "chained")
  methods_to_run <- if (identical(type, "all")) all_methods else intersect(type, all_methods)

  if (length(methods_to_run) == 0) {
    stop("No valid methods specified in 'type'.")
  }

  # Initialize lists to store results
  results_list <- list()
  equated_scores_list <- list()
  raw_scores_x <- seq(from = min_x, to = max_x, by = inc_x)

  # --- Tucker Method ---
  if ("tucker" %in% methods_to_run) {
    gamma1_tuck <- covxv1 / (sdv1^2)
    gamma2_tuck <- covyv2 / (sdv2^2)
    res_tuck <- linear_equate_observed_ci(mnx1, sdx1, mnv1, sdv1, mny2, sdy2, mnv2, sdv2, w1, mean_only, gamma1_tuck, gamma2_tuck)
    results_list$Tucker <- list(a = res_tuck$a, b = res_tuck$b)
    equated_scores_list$Tucker <- res_tuck$a * raw_scores_x + res_tuck$b
  }

  # --- Levine Observed/True Score Method ---
  # These are linked because they share gamma calculations
  if (any(c("levine_observed", "levine_true") %in% methods_to_run)) {
    if (anchor != 0) { # Internal anchor
      gamma1_lev <- sdx1^2 / covxv1
      gamma2_lev <- sdy2^2 / covyv2
    } else { # External anchor
      gamma1_lev <- (sdx1^2 + covxv1) / (sdv1^2 + covxv1)
      gamma2_lev <- (sdy2^2 + covyv2) / (sdv2^2 + covyv2)
    }

    if ("levine_observed" %in% methods_to_run) {
      res_lev_obs <- linear_equate_observed_ci(mnx1, sdx1, mnv1, sdv1, mny2, sdy2, mnv2, sdv2, w1, mean_only, gamma1_lev, gamma2_lev)
      results_list$`Levine Observed` <- list(a = res_lev_obs$a, b = res_lev_obs$b)
      equated_scores_list$`Levine Observed` <- res_lev_obs$a * raw_scores_x + res_lev_obs$b
    }

    if ("levine_true" %in% methods_to_run) {
      a_lev_true <- if (!mean_only) gamma2_lev / gamma1_lev else 1.0
      b_lev_true <- (mny2 - a_lev_true * mnx1) + gamma2_lev * (mnv1 - mnv2)
      results_list$`Levine True` <- list(a = a_lev_true, b = b_lev_true)
      equated_scores_list$`Levine True` <- a_lev_true * raw_scores_x + b_lev_true
    }
  }

  # --- Chained Linear Method ---
  if ("chained" %in% methods_to_run) {
    gamma1_chain <- sdx1 / sdv1
    gamma2_chain <- sdy2 / sdv2
    res_chain <- linear_equate_observed_ci(mnx1, sdx1, mnv1, sdv1, mny2, sdy2, mnv2, sdv2, w1, mean_only, gamma1_chain, gamma2_chain)
    results_list$Chained <- list(a = res_chain$a, b = res_chain$b)
    equated_scores_list$Chained <- res_chain$a * raw_scores_x + res_chain$b
  }

  # --- Combine results ---
  summary_df <- do.call(rbind, lapply(names(results_list), function(name) {
    data.frame(Method = name, a = results_list[[name]]$a, b = results_list[[name]]$b)
  }))

  equated_df <- as.data.frame(do.call(rbind, equated_scores_list))
  colnames(equated_df) <- paste0("x_", raw_scores_x)
  rownames(equated_df) <- names(equated_scores_list)

  return(list(summary = summary_df, equated_scores = equated_df))
}
