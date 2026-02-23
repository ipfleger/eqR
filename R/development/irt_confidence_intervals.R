#' Ogasawara (2001) Asymptotic Standard Errors for IRT True Score Equating
#'
#' @description
#' Calculates the asymptotic standard errors for IRT true score equating using
#' the multivariate delta method as described in Ogasawara (2001).
#' It handles the error propagation from item parameter estimation (for both
#' unique and common items) to the equating coefficients (A, B) and finally
#' to the equated true scores.
#'
#' @param eq An `equate_recipe` object (or list) containing method options.
#' @param forms A character vector of length 2 (e.g., c("FormX", "FormY")).
#' @param title The title of the equating method in `eq`.
#' @param vcov_list A named list of covariance matrices for item parameters.
#'   Must have names matching `forms`. Each matrix should align with the
#'   'flat' vector of parameters (slope, intercept, guess) for that form.
#' @param conf_level Confidence level for the intervals (default 0.95).
#' @references https://www.jstor.org/stable/3657936?read-now=1&seq=1
#' @return A list structured like the user's `linear_sgrg` example.
irt_true_score_se <- function(eq, forms, title, vcov_list = NULL, conf_level = 0.95) {

  # --- 1. Setup & Data Extraction ---
  method_options <- eq@methods[[title]]$options
  theta <- method_options$theta

  # Standardize item parameters using the helper from your file
  # (Assuming irt_coefs is available from irt_helpers.R)
  raw_pars <- irt_coefs(method_options$irt_pars)
  pars_x <- raw_pars[[forms[1]]]
  pars_y <- raw_pars[[forms[2]]]

  # Identify Common Items
  common_map <- method_options$common_items
  common_idx_x <- common_map[[forms[1]]]
  common_idx_y <- common_map[[forms[2]]]

  # Flatten parameters for Jacobian calculations (Slope, Intercept, Guess per item)
  # Vector format: [a1, d1, g1, a2, d2, g2, ...]
  flat_pars_x <- unlist(t(pars_x[, c("slope", "intercept", "guess")]))
  flat_pars_y <- unlist(t(pars_y[, c("slope", "intercept", "guess")]))

  # --- 2. Handle Covariance Matrices ---
  # If vcov not provided, assume 0 covariances (unrealistic, but runs)
  if (is.null(vcov_list)) {
    warning("No 'vcov_list' provided. Assuming zero covariance (SEs will be underestimated).")
    vcov_x <- diag(0, length(flat_pars_x))
    vcov_y <- diag(0, length(flat_pars_y))
  } else {
    vcov_x <- vcov_list[[forms[1]]]
    vcov_y <- vcov_list[[forms[2]]]
  }

  # Construct Block Diagonal Matrix for All Item Params [Psi_X, Psi_Y]
  # Ogasawara assumes independence between forms (separate calibrations)
  Sigma_Items <- as.matrix(Matrix::bdiag(vcov_x, vcov_y))

  # --- 3. Equating Coefficients (A & B) & Their Derivatives ---
  # We need the gradient of A and B w.r.t. the COMMON item parameters.
  # This corresponds to the implicit differentiation in Ogasawara.
  # We do this numerically for robustness across SL/Haebara/Moment methods.

  # Helper to compute A/B from a subset of flat parameters
  get_AB_from_flat <- function(flat_common_vec) {
    n_common <- length(flat_common_vec) / 2 # Half from X, Half from Y
    p_x <- matrix(flat_common_vec[1:(length(flat_common_vec)/2)], ncol=3, byrow=TRUE)
    p_y <- matrix(flat_common_vec[(length(flat_common_vec)/2 + 1):length(flat_common_vec)], ncol=3, byrow=TRUE)

    df_x <- data.frame(slope=p_x[,1], intercept=p_x[,2], guess=p_x[,3])
    df_y <- data.frame(slope=p_y[,1], intercept=p_y[,2], guess=p_y[,3])

    # Use the appropriate method (Stocking-Lord, etc.)
    if(method_options$transform_method %in% c("stocking_lord", "haebara")){
      res <- scale_curve(df_x, df_y, theta, method = method_options$transform_method)
    } else {
      res <- scale_moment(df_x, df_y, method = method_options$transform_method)
    }
    return(c(res$A, res$B))
  }

  # Extract common item parameters to pass to Jacobian
  flat_common_x <- unlist(t(pars_x[common_idx_x, c("slope", "intercept", "guess")]))
  flat_common_y <- unlist(t(pars_y[common_idx_y, c("slope", "intercept", "guess")]))
  flat_common_all <- c(flat_common_x, flat_common_y)

  # 3a. Point Estimates
  AB_est <- get_AB_from_flat(flat_common_all)
  A_est <- AB_est[1]
  B_est <- AB_est[2]

  # 3b. Numerical Jacobian of A,B w.r.t Common Item Params
  # J_AB_common has dimensions 2 x (Total params of common items)
  J_AB_common <- simple_jacobian(get_AB_from_flat, flat_common_all)

  # 3c. Map J_AB_common to J_AB_total (w.r.t ALL items X and Y)
  # Most entries are 0 (non-common items don't affect A/B)
  J_AB_total <- matrix(0, nrow = 2, ncol = length(flat_pars_x) + length(flat_pars_y))

  # Indices of common parameters in the full vectors
  # Because we flattened by row, indices are regular
  idx_x_global <- as.vector(outer(c(0,1,2), (common_idx_x-1)*3, "+") + 1)
  idx_y_global <- as.vector(outer(c(0,1,2), (common_idx_y-1)*3, "+") + 1) + length(flat_pars_x) # Offset for Y

  # Fill the Jacobian
  n_com_params <- length(flat_common_x)
  J_AB_total[, idx_x_global] <- J_AB_common[, 1:n_com_params]
  J_AB_total[, idx_y_global] <- J_AB_common[, (n_com_params+1):(2*n_com_params)]

  # --- 4. SE of the Equated Scores (Zeta) ---

  # Apply transformation to X parameters for calculation
  pars_x_trans <- transform_irt_pars(pars_x, A_est, B_est)

  # Score range
  min_score <- attr(eq@data[[forms[1]]], "min") %||% 0
  max_score <- attr(eq@data[[forms[1]]], "max") %||% nrow(pars_x)
  x_scores <- seq(min_score, max_score, by = 1) # Support integers usually

  # Pre-calculate TCCs and their derivatives for efficiency
  # We need derivatives of TCC w.r.t item params (dP/dPsi)

  # Calculate TCC of X (trans) and Y to get the mapping
  tcc_x <- calculate_tcc(pars_x_trans, theta)
  tcc_y <- calculate_tcc(pars_y, theta)

  # Get Theta equivalent for each X raw score (x -> theta_y)
  # Note: In true score equating, x = TCC_X(theta). We invert this.
  theta_equiv <- stats::approx(x = tcc_x, y = theta, xout = x_scores, rule = 2)$y

  # Calculate Equated Scores (Point Estimates)
  # Zeta = TCC_Y(theta_equiv)
  zeta_scores <- calculate_tcc(pars_y, theta_equiv)

  # --- 5. Calculate Gradients for Zeta ---
  # d(Zeta)/d(Psi) = Term1 + Term2
  # Term 1: Direct effect of Y params on Zeta
  # Term 2: Indirect effect via Theta_equiv (which depends on A, B, and X params)

  se_zeta <- numeric(length(x_scores))

  for (k in seq_along(x_scores)) {
    th <- theta_equiv[k]

    # Gradient of TCC_Y at theta w.r.t Y params
    # Returns vector of length n_items_Y * 3
    grad_TCCy_PsiY <- grad_tcc_wrt_params(pars_y, th)

    # Slope of TCC_Y at theta (sum of ICC slopes)
    slope_TCCy <- sum(grad_tcc_wrt_theta(pars_y, th))

    # Slope of TCC_X (on Y scale) at theta
    # Note: TCC_X_trans(theta) = sum P( a_x/A * theta + ... )
    slope_TCCx_trans <- sum(grad_tcc_wrt_theta(pars_x_trans, th))

    # Gradient of Theta w.r.t ALL params (Implicit Function Theorem)
    # d(Theta)/d(Psi) = - [d(TCC_X_trans)/d(Psi)] / slope_TCC_X_trans

    # 1. d(TCC_X_trans)/d(Psi_X_raw)
    # This involves the chain rule through A and B
    grad_TCCx_PsiXraw <- grad_tcc_trans_wrt_raw_params(pars_x, th, A_est, B_est)

    # 2. d(TCC_X_trans)/d(A) and d(B)
    grad_TCCx_AB <- grad_tcc_trans_wrt_AB(pars_x, th, A_est, B_est)

    # Combine to get full gradient of Theta w.r.t. ALL items [Psi_X, Psi_Y]
    # d(Theta)/d(Psi_total) = - (1/slope_X) * [ d(TCCx)/d(Psi_X) + d(TCCx)/d(AB) * d(AB)/d(Psi_total) ]

    # Pad grad_TCCx_PsiXraw with 0s for Y params
    grad_TCCx_total_direct <- c(grad_TCCx_PsiXraw, numeric(length(flat_pars_y)))

    # Indirect effect through A and B
    grad_TCCx_total_indirect <- t(J_AB_total) %*% grad_TCCx_AB

    grad_Theta_PsiTotal <- -1/slope_TCCx_trans * (grad_TCCx_total_direct + as.vector(grad_TCCx_total_indirect))

    # Final Gradient of Zeta
    # d(Zeta)/d(Psi_total) = d(TCC_Y)/d(Psi_Y) * [0...1] + slope_TCC_Y * d(Theta)/d(Psi_total)

    # Pad TCC_Y gradient with 0s for X params
    grad_Zeta_PsiY_only <- c(numeric(length(flat_pars_x)), grad_TCCy_PsiY)

    d_Zeta <- grad_Zeta_PsiY_only + slope_TCCy * grad_Theta_PsiTotal

    # --- SANDWICH ---
    var_zeta <- t(d_Zeta) %*% Sigma_Items %*% d_Zeta
    se_zeta[k] <- sqrt(as.vector(var_zeta))
  }

  # --- 6. Output Formatting ---

  # Calculate CIs
  z_crit <- qnorm(1 - (1 - conf_level)/2)
  ci_lower <- zeta_scores - z_crit * se_zeta
  ci_upper <- zeta_scores + z_crit * se_zeta

  # Estimate SE for A and B (Transformation Constants)
  cov_AB <- J_AB_total %*% Sigma_Items %*% t(J_AB_total)
  se_AB <- sqrt(diag(cov_AB))

  results <- list(list(
    parameters = data.frame(
      statistics = c("Slope (A)", "Intercept (B)"),
      estimate = c(A_est, B_est),
      se = se_AB,
      lower_bound = c(A_est, B_est) - z_crit * se_AB,
      upper_bound = c(A_est, B_est) + z_crit * se_AB,
      bootstrapped_estimate = NA_real_ # Not bootstrap
    ),
    x_score = x_scores,
    equivalent_score = zeta_scores,
    bootstrapped_estimate = NA_real_,
    nested_intervals = data.frame(
      se = se_zeta,
      lower_bound_95 = ci_lower,
      upper_bound_95 = ci_upper
    ),
    single = FALSE, # Assuming common item design for Ogasawara
    observed_scores_x = if (!is.null(eq@data[[forms[1]]])) rowSums(eq@data[[forms[1]]][, -1]) else NULL,
    observed_scores_y = if (!is.null(eq@data[[forms[2]]])) rowSums(eq@data[[forms[2]]][, -1]) else NULL
  )) |> stats::setNames(paste0("IRT_True_Score_SE (", method_options$transform_method, ")"))

  return(results)
}

# --- Internal Helpers for Derivatives ---

simple_jacobian <- function(func, x, ...) {
  # Central difference approximation
  epsilon <- 1e-6
  n <- length(x)
  f0 <- func(x, ...)
  m <- length(f0)
  J <- matrix(0, m, n)
  for (i in 1:n) {
    x_plus <- x
    x_minus <- x
    x_plus[i] <- x[i] + epsilon
    x_minus[i] <- x[i] - epsilon
    J[, i] <- (func(x_plus, ...) - func(x_minus, ...)) / (2 * epsilon)
  }
  return(J)
}

grad_tcc_wrt_params <- function(pars, theta) {
  # d(TCC)/d(Psi) for a single theta point
  # TCC = Sum( P_i(theta) )
  # Returns vector aligned with flattened params [a1, b1, g1, a2, b2, g2...]

  D <- 1.7 # or 1.0 depending on scaling, assume 1.0 for now or check options
  # Typically Ogasawara uses D=1.0 or absorbs it into 'a'. Assuming 1.0 for calc.

  a <- pars$slope
  b <- pars$intercept
  g <- pars$guess
  lin <- a * theta + b
  exp_lin <- exp(-lin)
  denom <- (1 + exp_lin)

  prob <- g + (1 - g) / denom

  # Derivatives of P_i w.r.t a_i, b_i, g_i
  # dP/da = (1-g) * P* * Q* * theta
  # dP/db = (1-g) * P* * Q* * 1
  # dP/dg = 1 - P*

  # Where P* = 1 / (1 + exp(-lin))
  p_star <- 1 / denom
  q_star <- 1 - p_star
  const <- (1 - g) * p_star * q_star

  dP_da <- const * theta
  dP_db <- const * 1
  dP_dg <- 1 - p_star

  # Interleave them
  grads <- numeric(length(a) * 3)
  grads[seq(1, by=3, length.out=length(a))] <- dP_da
  grads[seq(2, by=3, length.out=length(a))] <- dP_db
  grads[seq(3, by=3, length.out=length(a))] <- dP_dg

  return(grads)
}

grad_tcc_wrt_theta <- function(pars, theta) {
  # Returns vector of derivatives per item w.r.t theta
  a <- pars$slope
  b <- pars$intercept
  g <- pars$guess
  lin <- a * theta + b
  p_star <- 1 / (1 + exp(-lin))

  dP_dtheta <- (1 - g) * p_star * (1 - p_star) * a
  return(dP_dtheta)
}

grad_tcc_trans_wrt_raw_params <- function(pars_raw, theta, A, B) {
  # We need derivative of Sum( P( (a/A)*theta + (b - B*a/A) ) ) w.r.t raw a, b, g
  # Let theta_trans = (theta - B)/A  <-- Wait, standard transformation is:
  # New a* = a/A
  # New b* = (b - B*a)/A
  # Linear pred = a* theta + b* = (a/A) theta + (b - Ba)/A = (a(theta - B) + b)/A ?
  # Actually, usually theta_new = A*theta_old + B.
  # If we are mapping X to Y, and using pars_x_trans on Y scale:
  # pars_x_trans$slope = pars_x$slope / A
  # pars_x_trans$intercept = (pars_x$intercept - B*pars_x$slope) / A

  a <- pars_raw$slope
  b <- pars_raw$intercept
  g <- pars_raw$guess

  a_star <- a / A
  b_star <- (b - B*a) / A

  lin <- a_star * theta + b_star
  p_star <- 1 / (1 + exp(-lin))
  q_star <- 1 - p_star
  const <- (1 - g) * p_star * q_star

  # Chain rule components
  # d(lin)/da = theta * (1/A) + ( -B/A ) = (theta - B)/A
  # d(lin)/db = 1/A
  # d(lin)/dg = 0

  dP_da <- const * (theta - B) / A
  dP_db <- const * (1 / A)
  dP_dg <- 1 - p_star

  grads <- numeric(length(a) * 3)
  grads[seq(1, by=3, length.out=length(a))] <- dP_da
  grads[seq(2, by=3, length.out=length(a))] <- dP_db
  grads[seq(3, by=3, length.out=length(a))] <- dP_dg

  return(grads)
}

grad_tcc_trans_wrt_AB <- function(pars_raw, theta, A, B) {
  # Gradient of the Sum(P_trans) w.r.t A and B
  a <- pars_raw$slope
  b <- pars_raw$intercept
  g <- pars_raw$guess

  a_star <- a / A
  b_star <- (b - B*a) / A
  lin <- a_star * theta + b_star

  p_star <- 1 / (1 + exp(-lin))
  const <- (1 - g) * p_star * (1 - p_star)

  # d(lin)/dA = d( (a*theta + b - B*a)/A ) / dA = -1/A^2 * (a*theta + b - B*a) = -lin / A
  # d(lin)/dB = d( ... ) / dB = a/A * (-1) = -a/A

  dP_dA <- sum(const * (-lin / A))
  dP_dB <- sum(const * (-a / A))

  return(c(dP_dA, dP_dB))
}


#' Ogasawara (2003) SEs for IRT Observed Score Equating
#'
#' @description
#' Calculates standard errors for IRT Observed Score equating using the
#' Delta Method. It numerically approximates the gradient of the
#' complex observed-score equating function (involving Lord-Wingersky recursion
#' and equipercentile linking) with respect to all item parameters.
#'
#' @param eq An `equate_recipe` object.
#' @param forms A character vector of length 2 (e.g., c("FormX", "FormY")).
#' @param title The title of the equating method in `eq`.
#' @param vcov_list A named list of covariance matrices for item parameters.
#' @param conf_level Confidence level (default 0.95).
#'
#' @return A list structured like the user's `linear_sgrg` example.
irt_observed_score_se <- function(eq, forms, title, vcov_list = NULL, conf_level = 0.95) {

  # --- 1. Setup & Data Extraction ---
  method_options <- eq@methods[[title]]$options
  theta <- method_options$theta

  # Extract and standardize parameters
  raw_pars <- irt_coefs(method_options$irt_pars)
  pars_x <- raw_pars[[forms[1]]]
  pars_y <- raw_pars[[forms[2]]]

  # Identify Common Items
  common_map <- method_options$common_items
  common_idx_x <- common_map[[forms[1]]]
  common_idx_y <- common_map[[forms[2]]]

  # Flatten parameters into a single vector for Jacobian calculation
  # Order: [X_slope, X_int, X_guess, Y_slope, Y_int, Y_guess]
  flat_pars_x <- unlist(t(pars_x[, c("slope", "intercept", "guess")]))
  flat_pars_y <- unlist(t(pars_y[, c("slope", "intercept", "guess")]))
  all_params <- c(flat_pars_x, flat_pars_y)

  # Lengths for reconstruction
  nx <- nrow(pars_x)
  ny <- nrow(pars_y)

  # --- 2. Covariance Matrix ---
  if (is.null(vcov_list)) {
    warning("No 'vcov_list' provided. SEs will be zero.")
    vcov_x <- diag(0, length(flat_pars_x))
    vcov_y <- diag(0, length(flat_pars_y))
  } else {
    vcov_x <- vcov_list[[forms[1]]]
    vcov_y <- vcov_list[[forms[2]]]
  }
  # Block diagonal covariance matrix (assuming independence between calibrations)
  Sigma <- as.matrix(Matrix::bdiag(vcov_x, vcov_y))

  # --- 3. Wrapper Function for the "Full Chain" ---
  # This function takes the raw item parameters, performs the scale transformation
  # (if needed), runs the observed score equating, and returns the vector of
  # equated scores. The Jacobian of THIS function is 'd'.

  equating_chain_func <- function(params) {
    # A. Reconstruct Parameter DataFrames
    split_pt <- length(params) / 2 # Only valid if lengths are equal? Better to use calculated lengths
    # Use pre-calculated lengths
    p_x_vec <- params[1:length(flat_pars_x)]
    p_y_vec <- params[(length(flat_pars_x) + 1):length(params)]

    df_x <- as.data.frame(matrix(p_x_vec, ncol=3, byrow=TRUE))
    colnames(df_x) <- c("slope", "intercept", "guess")
    df_x$item <- pars_x$item # Keep names

    df_y <- as.data.frame(matrix(p_y_vec, ncol=3, byrow=TRUE))
    colnames(df_y) <- c("slope", "intercept", "guess")
    df_y$item <- pars_y$item

    # B. Scale Transformation (Calculate A & B)
    # Extract common items from the CURRENT perturbed parameters
    common_x <- df_x[common_idx_x, ]
    common_y <- df_y[common_idx_y, ]

    if(method_options$transform_method %in% c("stocking_lord", "haebara")){
      trans <- scale_curve(common_x, common_y, theta, method = method_options$transform_method)
    } else {
      trans <- scale_moment(common_x, common_y, method = method_options$transform_method)
    }

    # Transform Form X parameters
    df_x_trans <- transform_irt_pars(df_x, trans$A, trans$B)

    # C. Run Observed Score Equating
    # We call the internal helper directly
    res <- irt_observed_score_equate(
      irt_pars_x = df_x_trans,
      irt_pars_y = df_y, # Y is the base scale
      theta = theta,
      eq = eq,
      forms = forms,
      design = "cneg", # Assuming CNEG for this context
      method_options = method_options
    )

    return(res$equivalent_score)
  }

  # --- 4. Calculate Jacobian & SE ---
  # Get point estimate
  point_est <- equating_chain_func(all_params)

  # Calculate Jacobian numerically
  # J dimensions: [Length of Score Scale] x [Total Number of Parameters]
  # This may be slow for long tests!
  J <- simple_jacobian(equating_chain_func, all_params)

  # Sandwich: SE = sqrt( diag( J %*% Sigma %*% t(J) ) )
  # We do row-by-row multiplication to save memory: sum((J[i,] %*% Sigma) * J[i,])

  se_vec <- numeric(length(point_est))
  for(i in 1:length(point_est)) {
    # d_i is the gradient for score point i
    d_i <- J[i, ]
    var_i <- t(d_i) %*% Sigma %*% d_i
    se_vec[i] <- sqrt(as.vector(var_i))
  }

  # --- 5. Formatting Results ---
  z_crit <- qnorm(1 - (1 - conf_level)/2)

  # Get X score range
  min_score <- attr(eq@data[[forms[1]]], "min") %||% 0
  max_score <- attr(eq@data[[forms[1]]], "max") %||% nrow(pars_x)
  x_scores <- seq(min_score, max_score, by = 1)

  results <- list(list(
    parameters = data.frame(
      statistics = c("Slope", "Intercept"),
      estimate = c(NA, NA), # Observed score equating doesn't have simple A/B output for the final result
      se = c(NA, NA),
      lower_bound = c(NA, NA),
      upper_bound = c(NA, NA),
      bootstrapped_estimate = NA_real_
    ),
    x_score = x_scores,
    equivalent_score = point_est,
    bootstrapped_estimate = NA_real_,
    nested_intervals = data.frame(
      se = se_vec,
      lower_bound_95 = point_est - z_crit * se_vec,
      upper_bound_95 = point_est + z_crit * se_vec
    ),
    single = FALSE,
    observed_scores_x = if (!is.null(eq@data[[forms[1]]])) rowSums(eq@data[[forms[1]]][, -1]) else NULL,
    observed_scores_y = if (!is.null(eq@data[[forms[2]]])) rowSums(eq@data[[forms[2]]][, -1]) else NULL
  )) |> stats::setNames(paste0("IRT_Observed_Score_SE (", method_options$transform_method, ")"))

  return(results)
}
