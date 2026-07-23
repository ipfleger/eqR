#' Standardize IRT Item Parameters
#'
#' @description
#' A helper function to extract and standardize IRT item parameters from various
#' sources into a consistent data frame format (`slope`, `intercept`, `guess`).
#'
#' @param irt_pars The item parameters, which can be a data.frame, a `mirt`
#'   model object, a `TAM` model object, or a list of these.
#'
#' @return A data frame or a list of data frames with standardized item parameters.
#' @keywords internal
irt_coefs <- function(irt_pars) {
  if (is.list(irt_pars) && !is.data.frame(irt_pars)) {
    return(lapply(irt_pars, irt_coefs))
  }

  if (is.data.frame(irt_pars)) {
    # Standardize from a generic data frame
    # Handles both a/b/c and slope/intercept/guess formats
    df <- data.frame(
      item = irt_pars$item,
      slope = irt_pars$slope %||% irt_pars$a,
      intercept = irt_pars$intercept %||% (-1 * (irt_pars$a %||% 0) * (irt_pars$b %||% 0)),
      guess = irt_pars$guess %||% irt_pars$c %||% 0
    )

  } else if ((attr(class(irt_pars), "package") %||% "not") == "mirt") {
    # Extract from a mirt model object
    if (!requireNamespace("mirt", quietly = TRUE)) {
      cli::cli_abort("The 'mirt' package is required to extract coefficients from a mirt object.")
    }
    coefs <- mirt::coef(irt_pars, simplify = TRUE)$items
    df <- data.frame(
      item = rownames(coefs),
      slope = coefs[, "a1"],
      intercept = coefs[, "d"],
      guess = coefs[, "g"]
    )

  } else if (any(grepl("tam", class(irt_pars)))) {
    # Extract from a TAM model object
    df <- irt_pars$item[, c("item", "B.Cat1.Dim1", "AXsi_.Cat1", "guess")]
    colnames(df) <- c("item", "slope", "intercept", "guess")

  } else {
    cli::cli_abort("The format of 'item_params' is not recognized. Please provide a data.frame, mirt object, TAM object, or a list of these.")
  }
  rownames(df) <- df$item
  return(df)
}


#' Perform IRT Scale Transformation using Moment Methods
#'
#' @description
#' Estimates scale transformation constants A (slope) and B (intercept) using
#' either the mean/sigma or mean/mean method. These methods match moments of
#' the item parameters of the common items.
#'
#' @details
#' The function uses the `slope`/`intercept` parameterization internally. The
#' classic `b` difficulty parameter is equivalent to `-intercept / slope`.
#' \itemize{
#'   \item **Mean/Sigma:** Matches the mean and standard deviation of the `b` parameters.
#'   \item **Mean/Mean:** Matches the mean of the `b` parameters and the mean of the `a` (slope) parameters.
#' }
#'
#' @param common_items_x A data frame of standardized IRT parameters (`slope`,
#'   `intercept`, `guess`) for the common items from the new form (Form X).
#' @param common_items_y A data frame of standardized IRT parameters for the
#'   common items from the old form (Form Y).
#' @param method Either `"mean_sigma"` or `"mean_mean"`.
#'
#' @return A list containing the estimated slope `A` and intercept `B`.
#' @keywords internal
scale_moment <- function(common_items_x, common_items_y, method = "mean_sigma") {
  # Calculate classic 'b' parameter for moment matching
  b_x <- -common_items_x$intercept / common_items_x$slope
  b_y <- -common_items_y$intercept / common_items_y$slope

  if (method == "mean_sigma") {
    A <- sd(b_y) / sd(b_x)
    B <- mean(b_y) - A * mean(b_x)
  } else if (method == "mean_mean") {
    A <- mean(common_items_x$slope) / mean(common_items_y$slope)
    B <- mean(b_y) - A * mean(b_x)
  } else {
    cli::cli_abort("Unknown moment method: '{method}'. Choose 'mean_sigma' or 'mean_mean'.")
  }

  return(list(A = A, B = B))
}


#' Perform IRT Scale Transformation using Characteristic Curve Methods
#'
#' @description
#' Estimates scale transformation constants A (slope) and B (intercept) using
#' optimization to minimize the difference between item or test characteristic
#' curves of the common items.
#'
#' @details
#' This function uses `stats::optim()` to find the `A` and `B` constants that
#' minimize a criterion function. The method determines the nature of this
#' function:
#' \itemize{
#'   \item **stocking_lord:** Minimizes the sum of squared differences between the
#'     **Test** Characteristic Curves (TCCs) of the common items.
#'   \item **haebara:** Minimizes the sum of squared differences between all
#'     **Item** Characteristic Curves (ICCs) of the common items.
#' }
#'
#' @param common_items_x Standardized IRT parameters for common items on Form X.
#' @param common_items_y Standardized IRT parameters for common items on Form Y.
#' @param theta A numeric vector for the `theta` grid.
#' @param method Either `"stocking_lord"` or `"haebara"`.
#'
#' @return A list containing the estimated slope `A` and intercept `B`.
#' @keywords internal
scale_curve <- function(common_items_x, common_items_y, theta, method = "stocking_lord") {
  criterion_function <- function(params) {
    A <- params[1]
    B <- params[2]

    transformed_x <- transform_irt_pars(common_items_x, A, B)

    if (method == "stocking_lord") {
      # TCC method: Compare the sum of the curves
      tcc_y <- calculate_tcc(common_items_y, theta)
      tcc_x_transformed <- calculate_tcc(transformed_x, theta)
      error <- sum((tcc_y - tcc_x_transformed)^2)
    } else if (method == "haebara") {
      # ICC method: Compare each item's curve individually
      icc_y <- calculate_icc(common_items_y, theta)
      icc_x_transformed <- calculate_icc(transformed_x, theta)
      error <- sum((icc_y - icc_x_transformed)^2)
    } else {
      cli::cli_abort("Unknown curve method: '{method}'. Choose 'stocking_lord' or 'haebara'.")
    }
    return(error)
  }

  initial_params <- c(1, 0)
  opt_result <- stats::optim(par = initial_params, fn = criterion_function, method = "BFGS")

  return(list(A = opt_result$par[1], B = opt_result$par[2]))
}

#' Apply Scale Transformation to IRT Parameters
#'
#' @description
#' Transforms a set of IRT parameters from one scale to another using the
#' provided A and B constants, based on the slope-intercept parameterization.
#'
#' @param irt_pars A data frame of standardized IRT parameters.
#' @param A The slope of the transformation.
#' @param B The intercept of the transformation.
#'
#' @return A data frame with the transformed parameters.
#' @keywords internal
transform_irt_pars <- function(irt_pars, A, B) {
  transformed_pars <- irt_pars
  # Transformation formulas for slope-intercept parameters
  transformed_pars$slope <- irt_pars$slope / A
  transformed_pars$intercept <- (irt_pars$intercept - (B * irt_pars$slope)) / A
  return(transformed_pars)
}


#' IRT Equating Wrapper
#'
#' @description
#' This function is the main internal dispatcher for IRT equating. It determines
#' the equating design, performs scale transformation if necessary (for CNEG),
#' and then calls the appropriate equating engine (`true_score` or `observed_score`).
#'
#' @inheritParams irt_true_score_equate
#' @param design A character string for the equating design ("sg" or "cneg").
#' @param type A character string specifying the IRT equating type.
#' @param eq The `equate_recipe` object.
#' @param title The unique title of the method being run.
#'
#' @keywords internal
irt <- function(forms, design, type, eq, title, ...) {
  # run_equating() passes design codes "S"/"R"/"CG"; the IRT internals branch on
  # the lowercase "cneg" for common-item nonequivalent groups. Normalize here so
  # the scale-transformation path is actually taken for the CG design.
  design <- if (toupper(design) %in% c("CG", "CNEG")) "cneg" else tolower(design)

  method_options <- eq@methods[[title]]$options
  theta <- method_options$theta
  all_irt_pars <- irt_coefs(method_options$irt_pars) # Standardize

  transformation_details <- NULL

  if (design == "cneg") {
    if (!is.list(all_irt_pars) || is.data.frame(all_irt_pars)) {
      cli::cli_abort("For CNEG IRT, 'irt_pars' must be a named list of parameter sets, one for each form.")
    }
    if (!all(forms %in% names(all_irt_pars))) {
      cli::cli_abort("The names in the 'irt_pars' list must match the form names: '{forms[1]}' and '{forms[2]}'.")
    }

    raw_pars_x <- all_irt_pars[[forms[1]]]
    raw_pars_y <- all_irt_pars[[forms[2]]]

    # common_items_map is now guaranteed to exist by the `equate` dispatcher
    common_items_map <- method_options$common_items

    common_x <- raw_pars_x[common_items_map[[forms[1]]], ]
    common_y <- raw_pars_y[common_items_map[[forms[2]]], ]

    scale_func <- if (method_options$transform_method %in% c("stocking_lord", "haebara")) {
      scale_curve
    } else {
      scale_moment
    }
    trans_consts <- scale_func(common_x, common_y, theta, method = method_options$transform_method)

    irt_pars_x <- transform_irt_pars(raw_pars_x, trans_consts$A, trans_consts$B)
    irt_pars_y <- raw_pars_y

    # --- Perform anchor test equating for diagnostics ---
    if(method_options$anchor_type == "true_score"){
      anchor_equating <- irt_true_score_equate(
        irt_pars_x = transform_irt_pars(common_x, trans_consts$A, trans_consts$B),
        irt_pars_y = common_y,
        theta = theta,
        eq = eq,
        forms = forms,
        min_score = 0,
        max_score = nrow(common_x)
      )
    } else {
      anchor_equating <- irt_observed_score_equate(
        irt_pars_x = transform_irt_pars(common_x, trans_consts$A, trans_consts$B),
        irt_pars_y = common_y,
        theta = theta,
        eq = eq,
        forms = forms,
        design = "cneg", # Use CNEG for anchor as well
        method_options = method_options,
        min_score = 0,
        max_score = nrow(common_x)
      )
    }

    transformation_details <- list(
      common_items = common_items_map,
      transform_method = method_options$transform_method,
      A = trans_consts$A,
      B = trans_consts$B,
      anchor_test_equating = anchor_equating
    )

  } else {
    if (is.data.frame(all_irt_pars)) {
      rownames(all_irt_pars) <- all_irt_pars$item
      irt_pars_x <- all_irt_pars[eq@forms[[forms[1]]], ]
      irt_pars_y <- all_irt_pars[eq@forms[[forms[2]]], ]
    } else if (is.list(all_irt_pars)) {
      irt_pars_x <- all_irt_pars[[forms[1]]]
      irt_pars_y <- all_irt_pars[[forms[2]]]
    } else {
      cli::cli_abort("Unsupported 'irt_pars' format for this design.")
    }
  }


  if (type == "true_score") {
    equating_result <- irt_true_score_equate(
      irt_pars_x = irt_pars_x,
      irt_pars_y = irt_pars_y,
      theta = theta,
      eq = eq,
      forms = forms
    )
    result_name <- "IRT True Score"
  } else if (type == "observed_score") {
    equating_result <- irt_observed_score_equate(
      irt_pars_x = irt_pars_x,
      irt_pars_y = irt_pars_y,
      theta = theta,
      eq = eq,
      forms = forms,
      design = design,
      method_options = method_options
    )
    result_name <- "IRT Observed Score"
  } else {
    cli::cli_abort("Unknown IRT equating type: '{type}'.")
  }

  # Add transformation details to diagnostics if they exist
  if (!is.null(transformation_details)) {
    equating_result$diagnostics$scale_transformation <- transformation_details
  }

  return(stats::setNames(list(equating_result), result_name))
}


#' Perform IRT True Score Equating
#'
#' @description
#' This function performs IRT true score equating, which links scores on two
#' test forms based on a common underlying ability (`theta`) scale.
#'
#' @details
#' This function first calculates the Test Characteristic Curve (TCC) for each
#' form. For each integer score on Form X, it finds the corresponding true score,
#' finds the theta level that produces that true score on Form X, and then finds
#' the true score on Form Y at that same theta level. For scores outside the
#' obtainable true score range, linear extrapolation is used.
#'
#' @param irt_pars_x A data frame of IRT item parameters for Form X (on the same scale as Y).
#' @param irt_pars_y A data frame of IRT item parameters for Form Y.
#' @param theta A numeric vector for the `theta` grid.
#' @param eq The `equate_recipe` object.
#' @param forms A character vector of the two forms being equated.
#' @param min_score Optional. The minimum score for the score scale. If NULL,
#'   it's derived from the `eq` object.
#' @param max_score Optional. The maximum score for the score scale.
#'
#' @return A list containing the equating table and diagnostic information. The
#'   diagnostics include the `tcc_table`.
#' @keywords internal
irt_true_score_equate <- function(irt_pars_x, irt_pars_y, theta, eq, forms, min_score = NULL, max_score = NULL) {
  # --- 1. Determine Score Range ---
  if (is.null(min_score) || is.null(max_score)) {
    min_score <- attr(eq@data[[forms[1]]], "min")
    max_score <- attr(eq@data[[forms[1]]], "max")
  }
  x_score <- min_score:max_score

  # --- 2. Calculate TCCs ---
  tcc_x <- calculate_tcc(irt_pars_x, theta)
  tcc_y <- calculate_tcc(irt_pars_y, theta)
  tcc_table <- data.frame(theta = theta, tcc_x = tcc_x, tcc_y = tcc_y)

  # --- 3. Create Raw-to-Raw Conversion Table ---
  theta_for_x <- stats::approx(x = tcc_x, y = theta, xout = x_score, rule = 2)$y
  equivalent_true_score <- stats::approx(x = theta, y = tcc_y, xout = theta_for_x, rule = 2)$y

  # --- 4. Extrapolate for Scores Outside True Score Range ---
  lb_x <- sum(irt_pars_x$guess)
  ub_x <- nrow(irt_pars_x)
  lb_y <- sum(irt_pars_y$guess)
  ub_y <- nrow(irt_pars_y)

  scores_below <- x_score < lb_x
  scores_above <- x_score > ub_x

  if (any(scores_below)) {
    if (lb_x > 0) { # Avoid division by zero
      equivalent_true_score[scores_below] <- lb_y / lb_x * x_score[scores_below]
    } else {
      equivalent_true_score[scores_below] <- lb_y
    }
  }
  if (any(scores_above)) {
    slope_above <- (ub_y - lb_y) / (ub_x - lb_x)
    equivalent_true_score[scores_above] <- slope_above * (x_score[scores_above] - ub_x) + ub_y
  }

  # --- 5. Assemble Standardized Output ---
  result <- list(
    x_score = x_score,
    equivalent_score = equivalent_true_score,
    observed_scores_x = if (!is.null(eq@data[[forms[1]]])) rowSums(eq@data[[forms[1]]][eq@forms[[forms[1]]]], na.rm = TRUE) else NULL,
    observed_scores_y = if (!is.null(eq@data[[forms[2]]])) rowSums(eq@data[[forms[2]]][eq@forms[[forms[2]]]], na.rm = TRUE) else NULL,
    diagnostics = list(tcc_table = tcc_table)
  )

  return(result)
}

#' Perform IRT Observed Score Equating
#'
#' @description
#' Performs IRT observed score equating for SG and CNEG designs.
#'
#' @details
#' This procedure derives the observed score distributions for each form from
#' the IRT model and a population ability distribution. For the CNEG design, it
#' creates a synthetic population based on the two groups. It then performs an
#' equipercentile equating on these two model-based distributions.
#'
#' @inheritParams irt_true_score_equate
#' @param design The equating design.
#' @param method_options A list of options for the method.
#' @keywords internal
irt_observed_score_equate <- function(irt_pars_x, irt_pars_y, theta, eq, forms, design, method_options, min_score = NULL, max_score = NULL) {
  # --- 1. Determine Score Range ---
  if (is.null(min_score) || is.null(max_score)) {
    min_score <- attr(eq@data[[forms[1]]], "min")
    max_score <- attr(eq@data[[forms[1]]], "max")
  }
  x_score <- min_score:max_score

  if (design == "cneg") {
    theta_dist_list <- method_options$theta_dist
    if (!is.list(theta_dist_list) || is.data.frame(theta_dist_list) || !all(forms %in% names(theta_dist_list))) {
      cli::cli_abort("For CNEG observed score IRT, 'theta_dist' must be a named list of distributions, one for each form.")
    }
    theta_dist_x_raw <- theta_dist_list[[forms[1]]] %||% "normal"
    theta_dist_y_raw <- theta_dist_list[[forms[2]]] %||% "normal"
    w1 <- method_options$w1

    theta_dist_x <- get_theta_dist(theta_dist_x_raw, theta)
    theta_dist_y <- get_theta_dist(theta_dist_y_raw, theta)

    cond_pmf_x_x <- sapply(theta, function(th) lord_wingersky_recursion(irt_pars_x, th))
    cond_pmf_y_y <- sapply(theta, function(th) lord_wingersky_recursion(irt_pars_y, th))

    f1x <- as.vector(cond_pmf_x_x %*% theta_dist_x)
    g2y <- as.vector(cond_pmf_y_y %*% theta_dist_y)

    f2x <- as.vector(cond_pmf_x_x %*% theta_dist_y)
    g1y <- as.vector(cond_pmf_y_y %*% theta_dist_x)

    marginal_pmf_x <- w1 * f1x + (1 - w1) * f2x
    marginal_pmf_y <- w1 * g1y + (1 - w1) * g2y

  } else { # Single Group Design
    theta_dist_raw <- method_options$theta_dist
    theta_dist <- get_theta_dist(theta_dist_raw, theta)

    cond_pmf_x <- sapply(theta, function(th) lord_wingersky_recursion(irt_pars_x, th))
    cond_pmf_y <- sapply(theta, function(th) lord_wingersky_recursion(irt_pars_y, th))

    marginal_pmf_x <- as.vector(cond_pmf_x %*% theta_dist)
    marginal_pmf_y <- as.vector(cond_pmf_y %*% theta_dist)
  }

  score_distributions <- data.frame(
    observed_score = x_score,
    pmf_x = marginal_pmf_x,
    pmf_y = marginal_pmf_y
  )
  names(score_distributions)[2:3] <- forms

  cum_pmf_x <- cumsum(marginal_pmf_x)
  cum_pmf_y <- cumsum(marginal_pmf_y)

  # Equate the two IRT model-based observed-score distributions with standard
  # equipercentile equating (perc_rank + perc_point), matching K&B and the rest
  # of the package -- rather than an ad hoc (possibly non-monotone) spline inverse.
  ns    <- length(x_score)
  min_s <- x_score[1]
  inc_s <- if (ns > 1) x_score[2] - x_score[1] else 1
  prd_x <- perc_rank(x = x_score, min = min_s, max = x_score[ns], inc = inc_s, crfd = cum_pmf_x)
  equivalent_score <- EquiEquate(nsy = ns, miny = min_s, incy = inc_s,
                                 crfdy = cum_pmf_y, nsx = ns, prdx = prd_x)

  result <- list(
    x_score = x_score,
    equivalent_score = equivalent_score,
    observed_scores_x = if (!is.null(eq@data[[forms[1]]])) rowSums(eq@data[[forms[1]]][eq@forms[[forms[1]]]], na.rm = TRUE) else NULL,
    observed_scores_y = if (!is.null(eq@data[[forms[2]]])) rowSums(eq@data[[forms[2]]][eq@forms[[forms[2]]]], na.rm = TRUE) else NULL,
    diagnostics = list(
      score_distributions = score_distributions
    )
  )

  return(result)
}

#' Helper to get theta distribution vector
#' @keywords internal
get_theta_dist <- function(theta_dist_option, theta) {
  if (is.character(theta_dist_option) && length(theta_dist_option) == 1) {
    dist <- switch(tolower(theta_dist_option),
                   "normal" = dnorm(theta),
                   "uniform" = rep(1, length(theta)),
                   cli::cli_abort("Unknown 'theta_dist' string: '{theta_dist_option}'.")
    )
  } else if (is.numeric(theta_dist_option)) {
    if (length(theta_dist_option) != length(theta)) {
      cli::cli_abort("'theta_dist' vector must have the same length as the 'theta' grid.")
    }
    dist <- theta_dist_option
  } else {
    cli::cli_abort("'theta_dist' must be a numeric vector or a string.")
  }
  return(dist / sum(dist))
}


#' Calculate the Test Characteristic Curve (TCC)
#'
#' @description
#' A helper function to compute the TCC (expected true score) from a set of
#' dichotomous item parameters for a given set of ability values.
#'
#' @param irt_pars A data frame of item parameters (`slope`, `intercept`, `guess`).
#' @param theta A numeric vector of ability values.
#'
#' @return A numeric vector of expected true scores, one for each theta value.
#' @keywords internal
calculate_tcc <- function(irt_pars, theta) {
  icc_matrix <- calculate_icc(irt_pars, theta)
  tcc <- colSums(icc_matrix)
  return(tcc)
}

#' Calculate Item Characteristic Curves (ICCs)
#'
#' @description
#' A helper function to compute the ICCs from a set of dichotomous item
#' parameters for a given set of ability values.
#'
#' @param irt_pars A data frame of item parameters (`slope`, `intercept`, `guess`).
#' @param theta A numeric vector of ability values.
#'
#' @return A matrix of probabilities, with items in rows and theta points in columns.
#' @keywords internal
calculate_icc <- function(irt_pars, theta) {
  linear_predictor <- outer(irt_pars$slope, theta, "*") + irt_pars$intercept
  prob_matrix <- 1 / (1 + exp(-linear_predictor))
  icc_matrix <- irt_pars$guess + (1 - irt_pars$guess) * prob_matrix
  return(icc_matrix)
}


#' Calculate Conditional Score Probabilities via Lord-Wingersky Algorithm
#'
#' @description
#' This function implements the Lord-Wingersky recursive algorithm to efficiently
#' calculate the probability of achieving each possible raw score on a test,
#' conditional on a single, specific ability level (`theta_point`).
#'
#' @details
#' This is a direct R translation of the `recursion` function from `IRTst.c`.
#'
#' @param irt_pars A data frame of item parameters for a complete test.
#' @param theta_point A single numeric value for the ability level.
#'
#' @return A numeric vector where the `k`-th element is the probability of
#'   obtaining a total score of `k-1`.
#' @keywords internal
lord_wingersky_recursion <- function(irt_pars, theta_point) {
  p_values <- irt_pars$guess + (1 - irt_pars$guess) / (1 + exp(-(irt_pars$slope * theta_point + irt_pars$intercept)))
  q_values <- 1 - p_values
  num_items <- nrow(irt_pars)
  P <- numeric(num_items + 1)
  P[1] <- 1.0

  for (i in 1:num_items) {
    for (k in (i + 1):1) {
      prob_score_k_minus_1_prev <- if (k > 1) P[k - 1] else 0
      P[k] <- P[k] * q_values[i] + prob_score_k_minus_1_prev * p_values[i]
    }
  }
  return(P)
}

