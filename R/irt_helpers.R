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



#' IRT Equating Wrapper
#'
#' @description
#' This function serves as a wrapper for different types of IRT equating,
#' such as true score and observed score equating. It retrieves the necessary
#' parameters from the `equate_recipe` object and calls the appropriate
#' low-level IRT equating engine.
#'
#' @param forms A character vector of length two indicating the forms to be equated.
#' @param design A character string for the equating design (e.g., "cg").
#' @param type A character string specifying the IRT equating type
#'   (e.g., "true_score").
#' @param eq The `equate_recipe` object containing all forms and methods.
#' @param title The unique title of the method being run, used to retrieve
#'   the correct options from the `eq` object.
#' @param ... Additional arguments, for forward compatibility.
#'
#' @keywords internal
#'
irt <- function(forms, design, type, eq, title, ...) {
  # --- 1. Retrieve and Process IRT options from the recipe ---
  method_options <- eq@methods[[title]]$options

  # Standardize the IRT parameters (from mirt object, etc.) into a data frame
  all_irt_pars <- irt_coefs(method_options$irt_pars)
  theta <- method_options$theta

  # --- 2. Select the correct item parameter sets for the forms ---
  # This logic handles cases where irt_pars is a single data frame
  # (for all items) or a list of data frames (one for each form).
  if (is.data.frame(all_irt_pars)) {
    # Set rownames to allow for subsetting by item ID
    rownames(all_irt_pars) <- all_irt_pars$item

    # Subset the items belonging to each form
    irt_pars_x <- all_irt_pars[eq@forms[[forms[1]]], ]
    irt_pars_y <- all_irt_pars[eq@forms[[forms[2]]], ]

  } else if (is.list(all_irt_pars)) {
    # Select the correct data frame from the list by form name
    irt_pars_x <- all_irt_pars[[forms[1]]]
    irt_pars_y <- all_irt_pars[[forms[2]]]
  }

  # --- 3. Dispatch to the correct IRT engine based on type ---
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
    # For observed score equating, retrieve and process the theta distribution
    theta_dist_option <- method_options$theta_dist

    if (is.character(theta_dist_option) && length(theta_dist_option) == 1) {
      theta_dist <- switch(tolower(theta_dist_option),
                           "normal" = dnorm(theta),
                           "uniform" = rep(1, length(theta)),
                           cli::cli_abort("Unknown 'theta_dist' string: '{theta_dist_option}'. Choose 'normal', 'uniform', or provide a numeric vector.")
      )
    } else if (is.numeric(theta_dist_option)) {
      if (length(theta_dist_option) != length(theta)) {
        cli::cli_abort("'theta_dist' vector must have the same length as the 'theta' grid.")
      }
      theta_dist <- theta_dist_option
    } else {
      cli::cli_abort("'theta_dist' must be a numeric vector or one of the strings 'normal' or 'uniform'.")
    }

    # Normalize the final distribution to ensure it sums to 1
    theta_dist <- theta_dist / sum(theta_dist)

    equating_result <- irt_observed_score_equate(
      irt_pars_x = irt_pars_x,
      irt_pars_y = irt_pars_y,
      theta = theta,
      theta_dist = theta_dist,
      eq = eq,
      forms = forms
    )
    result_name <- "IRT Observed Score"
  } else {
    cli::cli_abort("Unknown IRT equating type: '{type}'.")
  }

  # --- 4. Return results in a named list for consistency ---
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
#' form. It then creates a raw-score-to-raw-score conversion table by using
#' linear interpolation. For each integer score on Form X, it finds the
#' corresponding true score (TCC), and then finds the true score on Form Y that
#' corresponds to the same theta level. Finally, it interpolates to find the
#' non-integer raw score on Form Y that corresponds to that true score.
#'
#' @param irt_pars_x A data frame of IRT item parameters for Form X.
#' @param irt_pars_y A data frame of IRT item parameters for Form Y.
#' @param theta A numeric vector for the `theta` grid.
#' @param eq The `equate_recipe` object.
#' @param forms A character vector of the two forms being equated.
#'
#' @return A list with the standardized equating output structure.
#' @keywords internal
irt_true_score_equate <- function(irt_pars_x, irt_pars_y, theta, eq, forms) {
  # --- 1. Determine Score Range ---
  min_score <- attr(eq@data[[forms[1]]], "min")
  max_score <- attr(eq@data[[forms[1]]], "max")
  x_score <- min_score:max_score

  # --- 2. Calculate TCCs ---
  tcc_x <- calculate_tcc(irt_pars_x, theta)
  tcc_y <- calculate_tcc(irt_pars_y, theta)
  tcc_table <- data.frame(theta = theta, tcc_x = tcc_x, tcc_y = tcc_y)

  # --- 3. Create Raw-to-Raw Conversion Table ---
  # Interpolate to find the equivalent score on Form Y for each integer score on X
  equivalent_score <- stats::approx(x = tcc_x, y = tcc_y, xout = x_score, rule = 2)$y


  # --- 4. Assemble Standardized Output ---
  result <- list(
    x_score = x_score,
    equivalent_score = equivalent_score,
    observed_scores_x = rowSums(eq@data[[forms[1]]][, -1]),
    observed_scores_y = rowSums(eq@data[[forms[2]]][, -1]),
    diagnostics = list(
      tcc_table = tcc_table
    )
  )

  return(result)
}

#' Perform IRT Observed Score Equating
#'
#' @description
#' This function performs IRT observed score equating using the method described
#' by Lord (1965).
#'
#' @details
#' The procedure derives the observed score distributions for each form from the
#' IRT model and a population ability distribution. It then performs an
#' equipercentile equating on these two model-based distributions.
#'
#' @param irt_pars_x A data frame of IRT item parameters for Form X.
#' @param irt_pars_y A data frame of IRT item parameters for Form Y.
#' @param theta A numeric vector for the `theta` grid.
#' @param theta_dist A numeric vector of weights for the `theta` grid.
#' @param eq The `equate_recipe` object.
#' @param forms A character vector of the two forms being equated.
#'
#' @return A list with the standardized equating output structure.
#' @keywords internal
irt_observed_score_equate <- function(irt_pars_x, irt_pars_y, theta, theta_dist,
                                      eq, forms) {
  # --- 1. Determine Score Range ---
  min_score <- attr(eq@data[[forms[1]]], "min")
  max_score <- attr(eq@data[[forms[1]]], "max")
  x_score <- min_score:max_score

  # --- 2. Calculate Conditional PMFs ---
  cond_pmf_x_matrix <- sapply(theta, function(th) {
    lord_wingersky_recursion(irt_pars = irt_pars_x, theta_point = th)
  })
  cond_pmf_y_matrix <- sapply(theta, function(th) {
    lord_wingersky_recursion(irt_pars = irt_pars_y, theta_point = th)
  })

  # --- 3. Integrate to get Marginal PMFs ---
  marginal_pmf_x <- as.vector(cond_pmf_x_matrix %*% theta_dist)
  marginal_pmf_y <- as.vector(cond_pmf_y_matrix %*% theta_dist)
  score_distributions <- data.frame(
    observed_score = x_score,
    pmf_x = marginal_pmf_x,
    pmf_y = marginal_pmf_y
  )
  names(score_distributions)[2:3] <- forms


  # --- 4. Perform Equipercentile Equating ---
  cum_pmf_x <- cumsum(marginal_pmf_x)
  cum_pmf_y <- cumsum(marginal_pmf_y)

  equivalent_score <- stats::spline(x = cum_pmf_y, y = x_score, xout = cum_pmf_x, method = "natural")$y

  # --- 5. Assemble Standardized Output ---
  result <- list(
    x_score = x_score,
    equivalent_score = equivalent_score,
    observed_scores_x = rowSums(eq@data[[forms[1]]][, -1]),
    observed_scores_y = rowSums(eq@data[[forms[2]]][, -1]),
    diagnostics = list(
      score_distributions = score_distributions
    )
  )

  return(result)
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
  # Vectorized calculation of the TCC
  linear_predictor <- outer(irt_pars$slope, theta, "*") + irt_pars$intercept
  prob_matrix <- 1 / (1 + exp(-linear_predictor))
  icc_matrix <- irt_pars$guess + (1 - irt_pars$guess) * prob_matrix
  tcc <- colSums(icc_matrix)
  return(tcc)
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
  # Calculate p-values (prob of correct response) for ALL items at the single theta point
  p_values <- irt_pars$guess + (1 - irt_pars$guess) / (1 + exp(-(irt_pars$slope * theta_point + irt_pars$intercept)))
  q_values <- 1 - p_values

  # Lord-Wingersky recursion
  num_items <- nrow(irt_pars)
  # P[k] will store the probability of getting a score of k-1
  P <- numeric(num_items + 1)
  P[1] <- 1.0 # Prob of score 0 after 0 items is 1

  for (i in 1:num_items) {
    # Loop backwards from the current max possible score down to 0
    for (k in (i + 1):1) {
      prob_score_k_minus_1_prev <- if (k > 1) P[k - 1] else 0
      P[k] <- P[k] * q_values[i] + prob_score_k_minus_1_prev * p_values[i]
    }
  }
  return(P)
}

