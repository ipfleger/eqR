# Title: R Functions for Printing Equating Results
#
# Description:
# This script contains user-facing functions designed to print and summarize
# the results from various equating procedures in a clear and readable format,
# similar to the output from the original Equating Recipes C library.

# Load required packages
# install.packages("cli") # Run if not installed

#' Print a Summary of Continuized Log-Linear Equating for CG Design
#'
#' @description
#' This function takes the results from a continuized log-linear equating run
#' for a Common-Item Non-Equivalent Groups (CG) design and prints a formatted
#' summary to the console. It is an R equivalent of the C library's `Print_CC`
#' function.
#'
#' @param pdata A list object (like the one returned by `smooth_bll`) containing
#'   the setup and parameter information for the analysis. It should include
#'   details about the smoothing model (e.g., polynomial degrees).
#' @param conversion_table A list or data frame containing the equated scores. It is
#'   expected to have a column for the original raw scores and additional columns
#'   for the equated scores from each method (e.g., "FE", "ChainedE").
#' @param title A character string for the main title of the output.
#'
#' @return This function does not return a value. It prints its output to the
#'   console.
#'
#' @examples
#' # This is a conceptual example of how you would use it.
#' # Assume 'pdata_list' and 'results_df' are returned from your main wrapper.
#'
#' # print_cll_cg(
#' #   title = "My Equating Run: Form B to Form A",
#' #   pdata = pdata_list,
#' #   conversion_table = results_df
#' # )
print_cll_cg <- function(title, eq, pdata, conversion_table) {

  # --- 1. Header Section ---
  cli::cli_h1(title)
  cli::cli_text("CLL Equating with CINEG Design (w1 = {pdata$wts})")
  cli::cli_rule()

  # --- 2. Model Information Section ---
  cli::cli_h2("Log-Linear Model Specification")
  cli::cli_ul()
  cli::cli_li("Form X (from-form) polynomial degree: {pdata[[1]]$cu}")
  cli::cli_li("Form Y (to-form) polynomial degree: {pdata[[2]]$cu}")
  cli::cli_li("Anchor Test (V) polynomial degree (X group): {pdata[[1]]$cv}")
  cli::cli_li("Anchor Test (V) polynomial degree (Y group): {pdata[[2]]$cv}")
  cli::cli_li("Number of cross-product terms: {pdata[[1]]$cuv}")
  cli::cli_end()

  # --- 3. Equated Scores Table ---
  cli::cli_h2("Equated Raw Scores")



  # Use print() with options for nice formatting
  print(
    conversion_table,
    row.names = FALSE,
    digits = 4
  )

  # --- 4. Moments Table ---
  cli::cli_rule()
  cli::cli_text("Moments of Equated Score Distributions")


  raw_from_scores <- rowSums(eq@data[[names(pdata[1])]][!colnames(eq@data[[names(pdata[1])]]) %in% get_anchors(eq)[[paste0(names(pdata), collapse =";")]]], na.rm = TRUE)
  raw_to_scores <- rowSums(eq@data[[names(pdata[2])]][!colnames(eq@data[[names(pdata[2])]]) %in% get_anchors(eq)[[paste0(names(pdata), collapse =";")]]], na.rm = TRUE)
  convert <- function(conversion_table){
    function(x)   conversion_table[x+1,]
  }
  converted_scores = convert(conversion_table)(raw_from_scores)
  # Calculate moments for each column of equated scores
  # moments_list <- c(list(raw_to_scores =  get_moments(scores = raw_to_scores,
  #                                                     rel_freq = raw_to_scores/sum(raw_to_scores))),
  #                        lapply(converted_scores[,-grep("^se$|bound", colnames(converted_scores))], function(col) {
  #   get_moments(scores = converted_scores$score, rel_freq = col / sum(col, na.rm = TRUE))
  # }))

  moments_df <- do.call(rbind, list(raw_to_scores = raw_to_scores,
                                    raw_from_scores = raw_from_scores,
                                    converted_from_scores_bs=converted_scores$equivalent,
                                    converted_from_scores_no_bs=converted_scores$equivalent_no_bs) |> lapply(moments))#as.data.frame(do.call(rbind, moments_list))

  print(
    moments_df,
    digits = 4
  )

  plot_equivalent(conversion_table, relative = FALSE)

  cli::cli_rule()
}

moments <- function(x){
  c(mean = mean(x), sd = sd(x), skewness = moments::skewness(x), kurtosis = moments::kurtosis(x))
}

print_linear_sgrg <- function(title, pdata) {

  # --- 1. Header Section ---
  cli::cli_h1(title)
  cli::cli_text("Linear Equating with Single Group/Random Group Design")
  cli::cli_rule()

  # --- 2. Model Information Section ---
  cli::cli_h2("Linear Model Coefficients")
  cli::cli_ul()
  cli::cli_li("Intercept: {round(pdata$parameters$estimate[pdata$parameters$statistics == 'Intercept'], 2)} ({paste0(round(pdata$parameters[pdata$parameters$statistics == 'Intercept',c('lower_bound_95','upper_bound_95')], 2), collapse = ', ')})")
  cli::cli_li("Slope: {round(pdata$parameters$estimate[pdata$parameters$statistics == 'Slope'], 2)} ({paste0(round(pdata$parameters[pdata$parameters$statistics == 'Slope',c('lower_bound_95','upper_bound_95')], 2), collapse = ', ')})")
  cli::cli_end()

  # --- 3. Equated Scores Table ---
  cli::cli_h2("Equated Raw Scores")


conversion_table <- data.frame(score = pdata$x_score, equivalent_score = pdata$equivalent_score,
                               bootstrapped_estimate = pdata$bootstrapped_estimate, pdata$nested_intervals)
  # Use print() with options for nice formatting
  print(
    conversion_table,
    row.names = FALSE,
    digits = 4
  )

  # --- 4. Moments Table ---
  cli::cli_rule()
  cli::cli_text("Moments of Equated Score Distributions")


  raw_from_scores <- pdata$observed_scores_x
  raw_to_scores <- pdata$observed_scores_y
  convert <- function(conversion_table){
    function(x)   conversion_table[x+1,]
  }
  converted_scores = convert(conversion_table)(raw_from_scores)
  # Calculate moments for each column of equated scores
  # moments_list <- c(list(raw_to_scores =  get_moments(scores = raw_to_scores,
  #                                                     rel_freq = raw_to_scores/sum(raw_to_scores))),
  #                        lapply(converted_scores[,-grep("^se$|bound", colnames(converted_scores))], function(col) {
  #   get_moments(scores = converted_scores$score, rel_freq = col / sum(col, na.rm = TRUE))
  # }))

  moments_df <- do.call(rbind, list(observed_y_scores = raw_to_scores,
                                    observed_x_scores = raw_from_scores,
                                    converted_x_scores = converted_scores$equivalent_score,
                                    converted_x_scores_boot = converted_scores$bootstrapped_estimate) |> lapply(moments))#as.data.frame(do.call(rbind, moments_list))

  print(
    moments_df,
    digits = 4
  )
  point_data <- data.frame(x = raw_from_scores, y = raw_to_scores)
  plot_equivalent(conversion_table, relative = FALSE, point_data = point_data)
  plot_equivalent(conversion_table, relative = TRUE, point_data = point_data)

  cli::cli_rule()

  return(c(pdata, list(moments = moments_df,
                       equate = convert,
                       converted_scores = converted_scores,
                       plots = list(score_conversion = \(gg = FALSE, relative = FALSE) plot_equivalent(conversion_table, relative = relative, point_data = point_data, gg = gg)))))
}

print_linear_cg <- function(title, pdata) {

  # --- 1. Header Section ---
  cli::cli_h1(title)
  cli::cli_text("Linear Equating with Common-Item Nonequivalent Groups Design")
  cli::cli_rule()

  # --- 2. Model Information Section ---
  cli::cli_h2("Linear Model Coefficients")
  cli::cli_ul()
  cli::cli_li("Intercept: {round(pdata$parameters$estimate[pdata$parameters$statistics == 'Intercept'], 2)} (95% CI: {paste0(round(pdata$parameters[pdata$parameters$statistics == 'Intercept', c('lower_bound_95','upper_bound_95')], 2), collapse = ', ')})")
  cli::cli_li("Slope: {round(pdata$parameters$estimate[pdata$parameters$statistics == 'Slope'], 2)} (95% CI: {paste0(round(pdata$parameters[pdata$parameters$statistics == 'Slope', c('lower_bound_95','upper_bound_95')], 2), collapse = ', ')})")
  cli::cli_end()

  # --- 3. Equated Scores Table ---
  cli::cli_h2("Equated Raw Scores")

  conversion_table <- data.frame(
    score = pdata$x_score,
    equivalent_score = pdata$equivalent_score,
    bootstrapped_estimate = pdata$bootstrapped_estimate,
    pdata$nested_intervals
  )

  # Use print() with options for nice formatting
  print(
    conversion_table,
    row.names = FALSE,
    digits = 4
  )

  # --- 4. Moments Table ---
  cli::cli_rule()
  cli::cli_text("Moments of Score Distributions")

  # Create a function to apply the equating transformation
  # This uses approxfun for linear interpolation between integer scores
  equate_fun <- approxfun(conversion_table$score, conversion_table$equivalent_score, rule = 2)

  # Apply the equating to the observed scores from the 'from' group
  converted_scores <- equate_fun(pdata$observed_scores_x)

  # Calculate moments for the original and equated distributions
  moments_df <- do.call(rbind, list(
    observed_y_scores = moments(pdata$observed_scores_y),
    observed_x_scores = moments(pdata$observed_scores_x),
    converted_x_scores = moments(converted_scores)
  ))

  print(
    moments_df,
    digits = 4
  )

  # --- 5. Plotting ---
  # # Assumes a 'plot_equivalent' function exists, similar to your sgrg version
  plot_equivalent(conversion_table, relative = FALSE, x_rug = pdata$observed_scores_x, y_rug = pdata$observed_scores_y)
  plot_equivalent(conversion_table, relative = TRUE, x_rug = pdata$observed_scores_x)

  cli::cli_rule()

  # --- 6. Return Enhanced Object ---
  # Return an object that includes the calculated moments and plotting functions
  c(pdata, list(
    moments = moments_df,
    equate = equate_fun,
    converted_scores = converted_scores,
    plots = list(
      score_conversion = \(gg = FALSE, relative = FALSE) plot_equivalent(conversion_table, relative = relative, gg = gg, x_rug = pdata$observed_scores_x, y = `if`(relative, NULL, pdata$observed_scores_y))
    )
  ))
}


#' Create a Summary Table for CINEG Linear Equating Methods
#'
#' @description
#' This function takes a list of results from the `linear_cg` function and
#' generates a formatted summary table comparing the estimated slope and
#' intercept for each equating method (Tucker, Levine, etc.).
#'
#' @param results_list A named list where each element is the result object
#'   for a single CINEG linear equating method.
#' @param title A character string for the main title of the summary output.
#'
#' @return A data frame containing the comparative summary statistics, which is
#'   also printed to the console in a formatted way.
#' @author Google's Gemini

summary.linear_cg <- function(results_list, title = "Comparative Summary of Linear Equating Methods") {

  # --- 1. Header Section ---
  cli::cli_h1(title)
  cli::cli_text("Comparison of estimated slope and intercept parameters across all methods.")
  cli::cli_rule()

  # --- 2. Extract and Combine Parameters ---
  # Use lapply to iterate through the list of results, extracting the parameters for each method.
  summary_df <- do.call(rbind, lapply(names(results_list), function(method_name) {
    # Extract the parameters data frame for the current method
    params <- results_list[[method_name]]$parameters

    # Get the slope and intercept estimates
    slope_est <- params$estimate[params$statistics == 'Slope']
    int_est <- params$estimate[params$statistics == 'Intercept']

    # Create a one-row data frame for this method
    data.frame(
      Method = method_name,
      Intercept = int_est,
      Slope = slope_est,
      stringsAsFactors = FALSE
    )
  }))

  # --- 3. Print the Formatted Table ---
  print(summary_df, row.names = FALSE, digits = 4)
  cli::cli_rule()

  # Return the summary data frame invisibly
  summary_df
}

