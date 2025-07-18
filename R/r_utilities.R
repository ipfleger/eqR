get_anchors <- function(eq){
  lapply(1:nrow(eq@plan) |> `names<-`(apply(eq@plan, 1, paste0,collapse = ";")), \(i) eq@forms[[eq@plan[i,"from"]]][eq@forms[[eq@plan[i,"from"]]] %in% eq@forms[[eq@plan[i,"to"]]]])
}

define_range <- function(eq, form, method_options){
  if (is.null(method_options$min_x)) {
    points <- apply(eq@data[[form]], 2, max, na.rm = TRUE)
    minx <- sum(apply(eq@data[[form]][eq@forms[[form]]], 2, min, na.rm = TRUE))
    inc <- median(points)
    maxx <- sum(points)
  } else {
    minx <- method_options$min_x
    maxx <- method_options$max_x
    inc <- method_options$inc_x
  }

  list(min = minx, max = maxx, inc = inc, range = seq(minx, maxx, by = inc))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Define and Summarize a Score Range
#'
#' @description
#' Calculates the theoretical minimum and maximum scores, and the score increment
#' based on item-level data. This function is flexible, allowing the user to
#' either have these parameters calculated automatically or to provide them as
#' manual overrides. It can also optionally print a summary of the range and
#' the observed total scores.
#'
#' @details
#' The function determines the score range by summing the minimum and maximum
#' possible scores for each item (column) in the `ctabs` data frame. The increment
#' is calculated as the median of the maximum possible scores for the items.
#' The `%||%` operator is used to provide a default value if an argument is `NULL`,
#' allowing for easy overriding of calculated parameters.
#'
#' @param ctabs A data frame containing item-level response data, where each
#'   column represents an item.
#' @param min_score (Optional) A numeric value to override the calculated minimum score.
#'   If `NULL` (default), the minimum is calculated from the data.
#' @param max_score (Optional) A numeric value to override the calculated maximum score.
#'   If `NULL` (default), the maximum is calculated from the data.
#' @param inc (Optional) A numeric value to override the calculated score increment.
#'   If `NULL` (default), the increment is calculated from the data.
#' @param verbose A logical value. If `TRUE`, prints a summary of the calculated
#'   range and observed score distribution to the console. Defaults to `FALSE`.
#'
#' @return A list containing the following elements:
#'   \item{min}{The final minimum score (either user-provided or calculated).}
#'   \item{max}{The final maximum score.}
#'   \item{inc}{The final score increment.}
#'   \item{range}{A numeric vector representing the full sequence of scores.}
#'   \item{num_sum}{A quantile summary of the observed total scores.}
#' @export
#'
#' @examples
#' # Create a sample item data frame
#' sample_ctabs <- data.frame(
#'   item1 = c(0, 1, 1, 0),
#'   item2 = c(0, 0, 1, 1),
#'   item3 = c(0, 1, 0, 1)
#' )
equate <- function(forms, method, design, type, eq, title,
                   boot_type = "perc", boot_replications = 1000){
  if(method == "L"){
    linear(forms = forms, design = design, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if(method == "E"){
    equipercentile(forms = forms, design = design, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if(method == "IRT"){
    irt(forms = forms, design = design, type = type, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  }
}


#'
#' # Get the range automatically and print the summary
#' range_info <- get_range(sample_ctabs, verbose = TRUE)
#'
#' # Get the range with a manual override for the minimum score
#' range_override <- get_range(sample_ctabs, min_score = 0)

get_range <- function(ctabs, form_name = NULL, min_score = NULL, max_score = NULL, inc = NULL, verbose = TRUE){

  # Calculate range parameters automatically from the data
  points <- apply(ctabs, 2, max, na.rm = TRUE)

  # Use user-provided values if they exist, otherwise use calculated defaults
  min_val <- min_score %||% sum(apply(ctabs, 2, min, na.rm = TRUE))
  inc <- inc %||% median(points)
  max_val <- max_score %||% sum(points)

  score_range <- seq(min_val, max_val, by = inc)

  # Calculate summary of observed total scores
  total_scores <- rowSums(ctabs, na.rm = TRUE)
  num_sum <- quantile(total_scores)

  # If verbose is TRUE, print the summary to the console
  if (verbose) {
    # Using the cli package for formatted output
    cli::cli_h1("Data Check:")
    cli::cli_inform("{form_name%||%'Score Summary'}: min = {min_val}, max = {max_val}, inc = {inc}")
    print(num_sum)
    cli::cli_rule()
  }

  # Return the list of calculated values
  list(min = min_val, max = max_val, inc = inc, range = score_range, num_sum = num_sum, points = points)
}

#' Print a Summary of Equating Results Object
#'
#' @param x An object of class `equate_results`.
#' @param ... Additional arguments (not used).
#' @return Invisibly returns the input object `x`.
#' @export
print.equate_results <- function(x, ...) {
  cli::cli_h1("Equating Results")
  cli::cli_text("This is an `equate_results` object. It contains the results of one or more equating plans.")
  cli::cli_rule()

  plan_names <- names(x)
  cli::cli_text("Found {length(plan_names)} equating plan{?s}: {.val {plan_names}}")

  for (plan_name in plan_names) {
    method_spec_names <- names(x[[plan_name]])
    cli::cli_text("Plan '{.strong {plan_name}}' contains {length(method_spec_names)} method specification{?s}: {.val {method_spec_names}}")
  }

  cli::cli_rule()
  cli::cli_alert_info("To see a detailed comparison of all methods for each plan, run {.fn summarize_equating_results} on this object.")

  invisible(x)
}


#' Summarize Equating Results
#'
#' Provides a comprehensive summary for each equating plan within an analysis
#' performed by `run_equating`. It generates comparative tables and plots for
#' the methods applied within each plan.
#'
#' @param results_object An object of class `equate_results`, the output from `run_equating`.
#' @param ... Additional arguments (not currently used).
#'
#' @return A nested list where each element corresponds to an equating plan.
#'   Each plan's element is a list containing the summary components
#'   (tables and plots) for that specific plan.
#' @export
summarize_equating_results <- function(results_object, ...) {

  all_summaries <- list()

  # --- Iterate through each equating plan in the results object ---
  for (plan_name in names(results_object)) {
    cli::cli_h1("Summary for Plan: {plan_name}")

    plan_result <- results_object[[plan_name]]

    # --- 1. Smart Data Wrangling for the CURRENT plan ---
    # Create a flat list of all results for this plan, with clean, unique names.
    results_list_for_plan <- build_smart_named_list(plan_result)


    # --- 2. Generate Summary Components for the CURRENT plan ---
    param_table <- compare_parameters(results_list_for_plan)
    conv_table <- create_conversion_table(results_list_for_plan)
    moments_table <- compare_moments(results_list_for_plan)
    plot_conv_gg <- plot_conversions(conv_table, gg = TRUE)
    plot_conv_base <- function() plot_conversions(conv_table, gg = FALSE)
    plot_dist_gg <- plot_distributions(results_list_for_plan, gg = TRUE)
    plot_dist_base <- function() plot_distributions(results_list_for_plan, gg = FALSE)

    # --- 3. Print Summaries to Console for the CURRENT plan ---
    if (!is.null(param_table)) {
      cli::cli_h2("Parameter Comparison (Linear Methods)")
      print(param_table)
      cli::cli_rule()
    }
    cli::cli_h2("Conversion Table")
    print(conv_table)
    cli::cli_rule()
    cli::cli_h2("Moment Comparison")
    print(moments_table)
    cli::cli_rule(left = "End of Summary for {plan_name}")


    # --- 4. Store the summary for the current plan ---
    all_summaries[[plan_name]] <- list(
      parameter_comparison = param_table,
      conversion_table = conv_table,
      moment_comparison = moments_table,
      conversion_plot_gg = plot_conv_gg,
      distribution_plot_gg = plot_dist_gg,
      conversion_plot_base = plot_conv_base,
      distribution_plot_base = plot_dist_base
    )
  }

  invisible(all_summaries)
}


# --- Helper Functions for summarize_equating_results ---

#' @noRd
build_smart_named_list <- function(plan_result) {
  # This function flattens the results list for a single plan and gives
  # each result a clean, descriptive, and unique name.

  temp_list <- list()
  # Create an intermediate list that tracks the parent spec for each result
  for (method_spec_name in names(plan_result)) {
    method_spec_result <- plan_result[[method_spec_name]]
    for (sub_method_name in names(method_spec_result)) {
      temp_list <- c(temp_list, list(list(
        spec_name = method_spec_name,
        sub_name = sub_method_name,
        result = method_spec_result[[sub_method_name]]
      )))
    }
  }

  # Get all the base sub-method names and find which ones are duplicated
  all_sub_names <- sapply(temp_list, `[[`, "sub_name")
  name_counts <- table(all_sub_names)

  final_list <- list()
  final_names <- sapply(temp_list, function(item) {
    sub_name <- item$sub_name
    # If the name is unique, just use it.
    if (name_counts[sub_name] == 1) {
      return(sub_name)
    }

    # If the name is duplicated, create a descriptive suffix.
    spec_name <- item$spec_name
    suffix_parts <- c()

    # Add suffix for mean_only
    if (grepl("mean_only", spec_name, fixed = TRUE)) {
      suffix_parts <- c(suffix_parts, "Mean Only")
    }

    # Add suffix for smoothing (if not the default 'N')
    spec_parts <- strsplit(spec_name, " ")[[1]]
    smooth_code <- spec_parts[4]
    if (smooth_code != "N") {
      smoother_name <- switch(smooth_code,
                              B = "Beta-Binomial", L = "Log-Linear", S = "Cubic Spline",
                              K = "Kernel", Z = "CLL", smooth_code)
      suffix_parts <- c(suffix_parts, paste("Smoothed:", smoother_name))
    }

    # If there are suffixes, combine them; otherwise, something else is duplicating the name.
    # Fallback to the full spec name to guarantee uniqueness if our suffix rules don't catch it.
    if (length(suffix_parts) > 0) {
      paste0(sub_name, " (", paste(suffix_parts, collapse = ", "), ")")
    } else {
      paste(spec_name, sub_name, sep = ".")
    }
  })

  # Assign the newly created names to the list of results
  results_list_for_plan <- lapply(temp_list, `[[`, "result")
  names(results_list_for_plan) <- final_names
  return(results_list_for_plan)
}


#' @noRd
compare_parameters <- function(results_list) {
  # Filter for methods that have a 'parameters' element
  linear_results <- Filter(function(x) !is.null(x$parameters), results_list)

  if (length(linear_results) == 0) {
    return(NULL)
  }

  param_list <- lapply(names(linear_results), function(method_name) {
    params <- linear_results[[method_name]]$parameters
    data.frame(
      Method = method_name,
      Intercept = params$estimate[params$statistics == 'Intercept'],
      Slope = params$estimate[params$statistics == 'Slope']
    )
  })

  do.call(rbind, param_list)
}

#' @noRd
create_conversion_table <- function(results_list) {
  # Extract the conversion from each result
  conversion_dfs <- lapply(names(results_list), function(method_name) {
    res <- results_list[[method_name]]
    df <- data.frame(
      x_score = res$x_score,
      value = res$equivalents
    )
    names(df)[2] <- method_name
    df
  })

  # Base R equivalent of purrr::reduce(..., dplyr::full_join)
  if (length(conversion_dfs) == 0) return(data.frame(x_score = numeric()))

  merged_df <- conversion_dfs[[1]]
  if (length(conversion_dfs) > 1) {
    for (i in 2:length(conversion_dfs)) {
      merged_df <- merge(merged_df, conversion_dfs[[i]], by = "x_score", all = TRUE)
    }
  }
  return(merged_df)
}

#' @noRd
compare_moments <- function(results_list) {
  # Function to safely calculate moments
  calculate_moments <- function(x) {
    if (length(x) < 2 || all(is.na(x))) return(c(mean=NA, sd=NA, skewness=NA, kurtosis=NA))
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      skewness = moments::skewness(x, na.rm = TRUE),
      kurtosis = moments::kurtosis(x, na.rm = TRUE)
    )
  }

  # Get the observed scores from the first result object (they are the same for all)
  obs_x <- results_list[[1]]$observed_scores_x
  obs_y <- results_list[[1]]$observed_scores_y

  # Calculate moments for observed scores
  moment_list <- list(
    Observed_From_Form = calculate_moments(obs_x),
    Observed_To_Form = calculate_moments(obs_y)
  )

  # Calculate moments for each set of converted scores
  for (method_name in names(results_list)) {
    res <- results_list[[method_name]]
    # Create an interpolation function for this method
    equate_fun <- stats::approxfun(res$x_score, res$equivalents, rule = 2)
    # Apply equating to the observed 'from' scores
    converted_scores <- equate_fun(obs_x)
    moment_list[[method_name]] <- calculate_moments(converted_scores)
  }

  # Combine into a data frame
  round(do.call(rbind, moment_list), 4)
}


#' @noRd
plot_conversions <- function(conversion_table, gg = TRUE) {
  # Base R equivalent of tidyr::pivot_longer
  plot_data_list <- lapply(names(conversion_table)[-1], function(method_name) {
    data.frame(
      x_score = conversion_table$x_score,
      Method = method_name,
      equivalent_score = conversion_table[[method_name]]
    )
  })
  plot_data <- do.call(rbind, plot_data_list)

  title <- "Comparison of Equating Functions"
  x_label <- "From-Form Score (X)"
  y_label <- "To-Form Score (Y)"

  if (gg) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_score, y = equivalent_score, color = Method, linetype = Method)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
      ggplot2::scale_color_viridis_d(option = "D") +
      ggplot2::labs(title = title, x = x_label, y = y_label) +
      ggplot2::theme_minimal()
    return(p)
  } else {
    # Base R plot
    x_range <- range(conversion_table$x_score, na.rm = TRUE)
    y_range <- range(conversion_table[, -1], na.rm = TRUE)
    plot(NA, xlim = x_range, ylim = y_range, xlab = x_label, ylab = y_label, main = title)
    abline(a = 0, b = 1, col = "grey50", lty = 2)

    colors <- viridis::viridis(ncol(conversion_table) - 1, option = "D")
    for (i in 2:ncol(conversion_table)) {
      lines(conversion_table$x_score, conversion_table[[i]], col = colors[i-1], lty = i-1, lwd = 2)
    }
    legend("topleft", legend = names(conversion_table)[-1], col = colors, lty = 1:(ncol(conversion_table)-1), lwd = 2, bty = "n")
  }
}

#' @noRd
plot_distributions <- function(results_list, gg = TRUE) {
  obs_y <- results_list[[1]]$observed_scores_y
  obs_x <- results_list[[1]]$observed_scores_x

  # Create a list of all score vectors to plot
  score_vectors <- list(
    Observed_To_Form = obs_y
  )
  for (method_name in names(results_list)) {
    res <- results_list[[method_name]]
    equate_fun <- stats::approxfun(res$x_score, res$equivalents, rule = 2)
    score_vectors[[method_name]] <- equate_fun(obs_x)
  }

  # Base R equivalent of purrr::map_df
  plot_data_list <- lapply(names(score_vectors), function(method_name) {
    data.frame(
      Method = method_name,
      score = score_vectors[[method_name]]
    )
  })
  plot_data <- do.call(rbind, plot_data_list)


  title <- "Comparison of Score Distributions"
  x_label <- "Score"

  if (gg) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = score, color = Method, linetype = Method)) +
      ggplot2::geom_density(linewidth = 1.1) +
      ggplot2::scale_color_viridis_d(option = "D") +
      ggplot2::labs(title = title, x = x_label) +
      ggplot2::theme_minimal()
    return(p)
  } else {
    # Base R plot
    densities <- lapply(score_vectors, density, na.rm = TRUE)
    x_range <- range(sapply(densities, function(d) d$x))
    y_range <- range(sapply(densities, function(d) d$y))

    plot(NA, xlim = x_range, ylim = y_range, xlab = x_label, ylab = "Density", main = title)
    colors <- viridis::viridis(length(densities), option = "D")
    for (i in seq_along(densities)) {
      lines(densities[[i]], col = colors[i], lty = i, lwd = 2)
    }
    legend("topright", legend = names(densities), col = colors, lty = 1:length(densities), lwd = 2, bty = "n")
  }
}
