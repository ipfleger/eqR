plot_equivalent <- function(data = NULL, results = NULL, relative = FALSE, gg = FALSE, point_data = NULL, x_rug = NULL, y_rug = NULL) {

  if(is.null(data)){
    if(!is.null(results)){
      data <- data.frame(score = results$x_score, equivalent_score = results$equivalent_score,
                                     bootstrapped_estimate = results$bootstrapped_estimate, results$nested_intervals)
      x_rug <- x_rug %||% results$observed_scores_x
      y_rug <- y_rug %||% results$observed_scores_y
    } else {
      cli::cli_abort("Must provide 'data' or 'results' to plot_equivalent")
    }
  }

  # Ensure the data is sorted by the x-axis variable
  data <- data[order(data$score), ]

  if(relative){
    # Define columns to modify
    cols_to_modify <- setdiff(names(data), "score")

    # Use vectorized subtraction
    data[cols_to_modify] <- data[cols_to_modify] - data$score

    if(!is.null(point_data)){
      point_data[[2]] <- point_data[[2]] - point_data[[1]]
    }


  }

  if(gg){
    p <- ggplot2::ggplot(data, ggplot2::aes(x = score, y= equivalent_score)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_bound_95, ymax = upper_bound_95), fill = "grey75", alpha = .5) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_bound_50, ymax = upper_bound_50), fill = "grey30", alpha = .5) +
      ggplot2::geom_abline(color = "firebrick", lty = 3, slope = !relative/1) +
      ggplot2::geom_line(ggplot2::aes(y = bootstrapped_estimate), color = "navy", lty = 2) +
      ggplot2::geom_line() +
      ggplot2::scale_y_continuous(breaks = scales::breaks_extended(20)) +
      ggplot2::scale_x_continuous(breaks = scales::breaks_extended(30))

    return(p)
  }
  # Define the full range for the x and y axes
  xlims <- range(data$score, na.rm = TRUE)
  ylims <- `if`(all(is.na(data$lower_bound_95)),
                range(data$equivalent_score, na.rm = TRUE),
                range(c(data$lower_bound_95, data$upper_bound_95), na.rm = TRUE))


  # 1. Set up the empty plot area
  plot(equivalent_score ~ score, data = data,
       type = "n",
       xlim = xlims,
       ylim = ylims,
       xlab = "score",
       ylab = "equivalent",
       main = ifelse(relative, "Relative Score Conversion Plot", "Score Conversion Plot"))

  # 2. Add light, non-distracting grid lines ✨
  grid(col = "gray85", lty = "solid", lwd = 0.5)

  if(!all(is.na(data$lower_bound_95))){
    # 3. Add the confidence ribbons
    polygon(x = c(data$score, rev(data$score)),
            y = c(data$lower_bound_95, rev(data$upper_bound_95)),
            col = rgb(0.75, 0.75, 0.75, alpha = 0.5),
            border = NA)

    polygon(x = c(data$score, rev(data$score)),
            y = c(data$lower_bound_50, rev(data$upper_bound_50)),
            col = rgb(0.30, 0.30, 0.30, alpha = 0.5),
            border = NA)
    lines(x = data$score, y = data$bootstrapped_estimate, col = "navy", lty = 2)
  }


  # 4. Add the various lines
  abline(a = 0, b = !relative/1, col = "firebrick", lty = 3)
  lines(x = data$score, y = data$equivalent_score)

  if(!is.null(point_data)){
    points(x = point_data[[1]], y = point_data[[2]],
           pch = 16, # Solid circle
           col = rgb(0, 0, 0, alpha = 0.2), # Semi-transparent black
           cex = 0.75)
  }

  # 5. Add rug plots if data is provided 🐾
  if (!is.null(x_rug)) {
    rug(x = x_rug, side = 1, col = rgb(0, 0, 0, 0.2), lwd = 1.5)
  }
  if (!is.null(y_rug)) {
    rug(x = y_rug, side = 2, col = rgb(0, 0, 0, 0.2), lwd = 1.5)
  }

}


#' Plot a Comparison of Multiple Equating Functions
#'
#' @description
#' This function takes a list of equating results and overlays their equating
#' functions on a single plot for easy visual comparison. It can generate
#' either a standard plot or a plot of the difference from the identity line.
#'
#' @param results_list A named list where each element is the result object
#'   for a single equating method. Each element must contain `x_score` and
#'   `equivalents`.
#' @param title A character string for the main title of the plot.
#' @param relative A logical value. If `TRUE`, the plot shows the difference
#'   between the equated score and the identity line (y - x). If `FALSE`
#'   (the default), it shows the standard equating function.
#' @param gg A logical value. If `TRUE`, the function returns a ggplot2 object.
#'   If `FALSE` (the default), it creates a base R plot.
#'
#' @return If `gg = FALSE`, the function is called for its side effect of creating a plot.
#'   If `gg = TRUE`, the function returns a ggplot object.
#' @author Google's Gemini

plot_equating_comparison <- function(results_list, title = NULL, relative = FALSE, gg = FALSE) {

  # --- Data Preparation ---
  # Combine all results into a single data frame for easier plotting
  plot_data_list <- lapply(names(results_list), function(method_name) {
    data.frame(
      Method = method_name,
      x_score = results_list[[method_name]]$x_score,
      equivalent_score = results_list[[method_name]]$equivalents
    )
  })
  plot_data <- do.call(rbind, plot_data_list)

  # If relative plot is requested, calculate the difference
  if (relative) {
    plot_data$equivalent_score <- plot_data$equivalent_score - plot_data$x_score
    y_label <- "Difference (Equated Y - X)"
    plot_title <- title %||% "Equating Function Differences from Identity"
    ref_line_intercept <- 0
  } else {
    y_label <- "Equated Form Y Score"
    plot_title <- title %||% "Comparison of Equating Functions"
    ref_line_intercept <- NA # Will be y=x instead
  }

  # --- ggplot2 version ---
  if (gg) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for ggplot output. Please install it.")
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_score, y = equivalent_score, color = Method, linetype = Method)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::labs(
        title = plot_title,
        x = "Form X Score",
        y = y_label,
        color = "Equating Method",
        linetype = "Equating Method"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    # Add the appropriate reference line
    if (relative) {
      p <- p + ggplot2::geom_hline(yintercept = ref_line_intercept, color = "gray70", linetype = "dashed")
    } else {
      p <- p + ggplot2::geom_abline(intercept = 0, slope = 1, color = "gray70", linetype = "dashed")
    }

    print(p)
    return(invisible(p))
  }

  # --- Base R version ---

  # 1. Setup Plotting Parameters
  colors <- c("black", "red", "blue", "green4")
  line_types <- 1:4

  x_range <- range(plot_data$x_score, na.rm = TRUE)
  y_range <- range(plot_data$equivalent_score, na.rm = TRUE)

  # 2. Initialize the Plot
  plot(
    NA, # Plot nothing initially
    xlim = x_range,
    ylim = y_range,
    xlab = "Form X Score",
    ylab = y_label,
    main = plot_title
  )

  # Add reference line
  if (relative) {
    abline(h = ref_line_intercept, col = "gray70", lty = 2)
  } else {
    abline(a = 0, b = 1, col = "gray70", lty = 2)
  }

  # 3. Overlay Functions
  for (i in 1:length(results_list)) {
    method_name <- names(results_list)[i]
    method_data <- plot_data[plot_data$Method == method_name, ]
    lines(
      method_data$x_score,
      method_data$equivalent_score,
      col = colors[i],
      lty = line_types[i],
      lwd = 2
    )
  }

  # 4. Add Legend
  legend(
    "topleft",
    legend = names(results_list),
    col = colors[1:length(results_list)],
    lty = line_types[1:length(results_list)],
    lwd = 2,
    bty = "n"
  )
}

# Helper to provide default titles
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
