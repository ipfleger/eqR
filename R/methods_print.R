#' Print an equate_recipe Object
#'
#' @description
#' Provides a summary of an `equate_recipe` object. The output adapts based on
#' whether the recipe has been run via `run_equating()`.
#'
#' - If the recipe has **not** been run, it prints a summary of the planned
#'   equatings, including the design and methods.
#' - If the recipe **has** been run, it prints a concise summary of the results
#'   for each equating plan and method.
#'
#' @param x An object of class `equate_recipe`.
#' @param ... Not used.
#'
#' @details Registered as an S7 method for the base `print()` generic via
#'   [S7::methods_register()] in `.onLoad`, so it needs no explicit NAMESPACE
#'   export.
#' @name print.equate_recipe
#' @keywords internal
S7::method(print, equate_recipe) <- function(x, ...) {
  if (length(x@results) == 0) {
    # --- Print the Plan (Recipe has NOT been run) ---
    cli::cli_h1("Equating Recipe")
    cli::cli_text("Design: {.val {x@design}}")

    cli::cli_h2("Forms")
    cli::cli_ul(items = names(x@forms))

    cli::cli_h2("Plan")
    for (i in 1:nrow(x@plan)) {
      cli::cli_text("{x@plan$from[i]} -> {x@plan$to[i]}")
    }

    cli::cli_h2("Methods")
    method_titles <- sapply(x@methods, function(m) translate_title(m$title))
    cli::cli_ul(items = method_titles)
  } else {
    # --- Print the Results (Recipe HAS been run) ---
    cli::cli_h1("Equating Results")
    cli::cli_text("Design: {.val {x@design}}")

    for (plan_name in names(x@results)) {
      plan_results <- x@results[[plan_name]]
      from_form <- strsplit(plan_name, ";")[[1]][1]
      to_form <- strsplit(plan_name, ";")[[1]][2]
      cli::cli_h2("Equating Plan: {from_form} to {to_form}")

      for (method_title in names(plan_results)) {
        cli::cli_rule(translate_title(method_title))
        result <- plan_results[[method_title]]
        print_result_summary(result, x, from_form)
      }
    }
    cli::cli_rule()
    cli::cli_alert_info("Tidy accessors: {.fn conversion_table}, {.fn conversions}, {.fn equated}; or {.code summary()} / {.code plot()}.")
  }
  invisible(x)
}


# --- Helper Functions for Printing ---

#' Print a concise summary of a single equating result
#' @keywords internal
print_result_summary <- function(result, recipe, from_form) {
  # Handle failed methods first
  if (all(is.na(result))) {
    cli::cli_alert_danger("Method failed to produce results.")
    return(invisible())
  }

  # The actual result might be nested one level down
  if (length(result) == 1 && is.list(result[[1]])) {
    result_name <- names(result)
    result <- result[[1]]
    attr(result, "result_name") <- result_name
  }


  # For IRT true score, the output is different
  if (is.data.frame(result) && all(c("theta", "tcc_x", "tcc_y") %in% names(result))) {
    cli::cli_text("IRT True Score Equating")
    cli::cli_ul(items = c(
      "TCC X Range: {round(min(result$tcc_x), 2)} to {round(max(result$tcc_x), 2)}",
      "TCC Y Range: {round(min(result$tcc_y), 2)} to {round(max(result$tcc_y), 2)}"
    ))
    return(invisible())
  }

  # For IRT observed score, use the conversion table
  if(is.list(result) && "conversion_table" %in% names(result)){
    result <- result$conversion_table
  }

  # For Linear methods, print parameters
  if (!is.null(result$parameters)) {
    cli::cli_text("Parameters:")
    print(knitr::kable(result$parameters, digits = 3))
  }

  # Summary of score conversions at key points
  obs_scores <- result$observed_scores_x
  if(is.null(obs_scores)) obs_scores <- recipe@data[[from_form]][[1]] # Fallback for IRT

  quantiles <- stats::quantile(obs_scores, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)

  # Ensure we have unique points to look up
  lookup_scores <- unique(round(quantiles))

  # Find the equated scores for these key points
  summary_df <- data.frame(
    Original_Score = lookup_scores,
    Equated_Score = result$equivalent_score[match(lookup_scores, result$x_score)]
  )

  if (!is.null(result$nested_intervals) && (is.matrix(result$nested_intervals) || is.data.frame(result$nested_intervals))) {
    summary_df$SE <- result$nested_intervals[match(lookup_scores, result$x_score), "se"]
  }


  cli::cli_text("\nScore Conversion Summary:")
  print(knitr::kable(summary_df, digits = 2))
  cat("\n")
}


#' Translate a method title into a human-readable string
#' @keywords internal
translate_title <- function(title) {
  parts <- strsplit(title, " ")[[1]]
  design <- switch(parts[1],
                   "S" = "Single Group",
                   "R" = "Random Groups",
                   "C" = "Common-Item",
                   parts[1]
  )
  method <- switch(parts[2],
                   "I" = "Identity",
                   "L" = "Linear",
                   "E" = "Equipercentile",
                   "IRT" = "IRT",
                   parts[2]
  )
  type <- tools::toTitleCase(gsub("_", " ", parts[3]))
  smooth <- switch(parts[4],
                   "N" = "(No Smoothing)",
                   "B" = "(Beta-Binomial)",
                   "L" = "(Log-Linear)",
                   "S" = "(Cubic Spline)",
                   "K" = "(Kernel)",
                   "Z" = "(CLL)",
                   ""
  )

  paste(design, ": ", method, " (", type, ") ", smooth, sep = "")
}


#' Summarize a run equating recipe
#'
#' For each form pairing, prints (and invisibly returns) the wide conversion
#' table and a moment comparison of the observed and equated distributions.
#'
#' @param object An `equate_recipe` that has been through [run_equating()].
#' @param ... Optional `from`/`to` to restrict to a single pairing.
#' @return Invisibly, a named list (one element per pairing) of
#'   `conversion_table` and `moments`.
#' @name summary.equate_recipe
#' @keywords internal
S7::method(summary, equate_recipe) <- function(object, ...) {
  conv <- object@conversions
  if (is.null(conv) || !nrow(conv)) {
    cli::cli_alert_warning("No results yet. Run {.fn run_equating} first.")
    return(invisible(NULL))
  }

  dots <- list(...)
  pairs <- unique(conv[, c("from", "to")])
  if (!is.null(dots$from)) pairs <- pairs[pairs$from == dots$from, , drop = FALSE]
  if (!is.null(dots$to))   pairs <- pairs[pairs$to   == dots$to,   , drop = FALSE]

  cli::cli_h1("Equating summary")
  out <- list()
  for (k in seq_len(nrow(pairs))) {
    from <- pairs$from[k]; to <- pairs$to[k]
    cli::cli_h2("{from} -> {to}")

    ct  <- conversion_table(object, from = from, to = to)
    mom <- .conversion_moments(object, from, to)

    cli::cli_text("Conversion table ({nrow(ct)} score points, {ncol(ct) - 1} method{?s}):")
    print(utils::head(ct, 10L), row.names = FALSE)
    if (nrow(ct) > 10L) cli::cli_text("... ({nrow(ct) - 10L} more rows; see {.fn conversion_table})")

    if (!is.null(mom)) {
      cli::cli_text("")
      cli::cli_text("Score-distribution moments:")
      print(round(mom, 3))
    }
    cli::cli_rule()
    out[[paste0(from, " -> ", to)]] <- list(conversion_table = ct, moments = mom)
  }
  invisible(out)
}


#' Plot equating conversions
#'
#' Draws the equated-score conversion line for each method of one form pairing,
#' with the identity line for reference (base graphics, no extra dependencies).
#'
#' @param x An `equate_recipe` that has been through [run_equating()].
#' @param ... Optional `from`/`to` to select the pairing.
#' @return Invisibly, the wide conversion table that was plotted.
#' @name plot.equate_recipe
#' @keywords internal
S7::method(plot, equate_recipe) <- function(x, ...) {
  dots <- list(...)
  ct <- conversion_table(x, from = dots$from, to = dots$to)
  methods <- setdiff(names(ct), "x_score")
  if (!length(methods)) cli::cli_abort("Nothing to plot.")

  cols <- grDevices::hcl.colors(max(length(methods), 2), palette = "viridis")
  xr <- range(ct$x_score, na.rm = TRUE)
  yr <- range(unlist(ct[methods]), na.rm = TRUE)

  graphics::plot(NA, xlim = xr, ylim = yr,
                 xlab = paste0(attr(ct, "from"), " raw score"),
                 ylab = paste0(attr(ct, "to"), " equivalent score"),
                 main = paste0("Equating: ", attr(ct, "from"), " -> ", attr(ct, "to")))
  graphics::abline(0, 1, col = "grey70", lty = 2)
  for (i in seq_along(methods)) {
    graphics::lines(ct$x_score, ct[[methods[i]]], col = cols[i], lwd = 2)
  }
  graphics::legend("topleft", legend = methods, col = cols[seq_along(methods)],
                   lwd = 2, bty = "n", cex = 0.9)
  invisible(ct)
}
