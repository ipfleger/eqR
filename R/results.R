# ------------------------------------------------------------------------------
# Tidy result layer for equate_recipe objects.
#
# The equating engines each return a nested list keyed by plan and method title,
# which is powerful but awkward to program against
# (`res@results$`Form C;Form A`$`S L all N`$Linear$equivalent_score`).
#
# run_equating() flattens all of that into a single tidy data frame stored in
# `@conversions`, and the accessors below are the primary way to use results:
#   - conversions()       long tidy table (one row per from/to/method/x_score)
#   - conversion_table()  wide, report-ready table (x_score + one col per method)
#   - equated()           apply a chosen conversion to arbitrary raw scores
# ------------------------------------------------------------------------------

#' Human-readable label for a smoothing code
#' @keywords internal
smooth_label <- function(code) {
  switch(code %||% "N",
         N = "none",  B = "beta-binomial",  L = "log-linear",
         S = "cubic spline",  K = "kernel",  Z = "continuized log-linear",
         code)
}

#' Flatten the raw nested results into one tidy conversion table
#'
#' Walks `@results` and returns a long data frame. Reads each method's stored
#' fields (never re-parses the space-delimited method title). Failed methods and
#' results whose shape has no `x_score`/`equivalent_score` (e.g. IRT true-score
#' curves) are skipped.
#' @keywords internal
collect_conversions <- function(eq) {
  if (!length(eq@results)) return(data.frame())

  parts <- list()
  for (plan_key in names(eq@results)) {
    frms <- strsplit(plan_key, ";", fixed = TRUE)[[1]]
    from <- frms[1]; to <- frms[2]
    plan_res <- eq@results[[plan_key]]
    if (!is.list(plan_res)) next

    for (title in names(plan_res)) {
      group <- plan_res[[title]]
      if (inherits(group, "equate_failed") || !is.list(group)) next
      smooth <- smooth_label(eq@methods[[title]]$smooth)

      for (submethod in names(group)) {
        r <- group[[submethod]]
        if (!is.list(r)) next
        xs <- r$x_score; es <- r$equivalent_score
        if (is.null(xs) || is.null(es) || length(xs) != length(es) || !length(xs)) next

        se <- rep(NA_real_, length(xs))
        ni <- r$nested_intervals
        if (!is.null(ni) && (is.data.frame(ni) || is.matrix(ni)) &&
            "se" %in% colnames(ni) && nrow(ni) == length(xs)) {
          se <- as.numeric(ni[, "se"])
        }

        parts[[length(parts) + 1L]] <- data.frame(
          from = from, to = to,
          method = submethod,
          smooth = smooth,
          x_score = as.numeric(xs),
          equivalent_score = as.numeric(es),
          se = se,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(parts)) return(data.frame())
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Resolve a single from/to pairing from a conversions table
#' @keywords internal
.pick_pair <- function(conv, from = NULL, to = NULL) {
  pairs <- unique(conv[, c("from", "to")])
  keep <- rep(TRUE, nrow(pairs))
  if (!is.null(from)) keep <- keep & pairs$from == from
  if (!is.null(to))   keep <- keep & pairs$to == to
  pairs <- pairs[keep, , drop = FALSE]
  if (!nrow(pairs)) {
    cli::cli_abort("No pairing matches from = {.val {from}}, to = {.val {to}}.")
  }
  if (nrow(pairs) > 1 && (is.null(from) || is.null(to))) {
    cli::cli_alert_info(
      "Multiple pairings available; using {.val {paste(pairs$from[1], '->', pairs$to[1])}}. Pass {.arg from}/{.arg to} to choose another."
    )
  }
  pairs[1, ]
}

#' Tidy equating results
#'
#' Returns the long, tidy conversion table produced by [run_equating()]: one row
#' per `from` form, `to` form, `method`, and raw score (`x_score`), with the
#' equated score and (when bootstrapped) its standard error.
#'
#' @param eq An `equate_recipe` that has been through [run_equating()].
#' @return A data frame with columns `from`, `to`, `method`, `smooth`,
#'   `x_score`, `equivalent_score`, and `se`.
#' @seealso [conversion_table()] for a wide, report-ready layout and [equated()]
#'   to convert arbitrary raw scores.
#' @export
conversions <- function(eq) {
  conv <- eq@conversions
  if (is.null(conv) || !nrow(conv)) {
    cli::cli_alert_warning("No conversions available. Did you run {.fn run_equating}?")
    return(invisible(data.frame()))
  }
  conv
}

#' Wide, report-ready conversion table
#'
#' Reshapes the tidy [conversions()] into a raw-score-by-method table for a
#' single form pairing: a `x_score` column followed by one column of equated
#' scores per method. This is the layout you typically paste into a report.
#'
#' @param eq An `equate_recipe` that has been through [run_equating()].
#' @param from,to Optional form names selecting the pairing. If omitted and only
#'   one pairing exists it is used; otherwise the first is used with a message.
#' @return A data frame with a `x_score` column and one column per method. The
#'   selected `from`/`to` are attached as attributes.
#' @export
conversion_table <- function(eq, from = NULL, to = NULL) {
  conv <- eq@conversions
  if (is.null(conv) || !nrow(conv)) {
    cli::cli_abort("No conversions available. Run {.fn run_equating} first.")
  }
  sel <- .pick_pair(conv, from, to)
  sub <- conv[conv$from == sel$from & conv$to == sel$to, , drop = FALSE]

  xs <- sort(unique(sub$x_score))
  wide <- data.frame(x_score = xs)
  for (m in unique(sub$method)) {
    mm <- sub[sub$method == m, , drop = FALSE]
    wide[[m]] <- mm$equivalent_score[match(xs, mm$x_score)]
  }
  rownames(wide) <- NULL
  attr(wide, "from") <- sel$from
  attr(wide, "to")   <- sel$to
  wide
}

#' Convert raw scores through a fitted equating
#'
#' Applies a chosen equating function to arbitrary raw scores on the `from`
#' form, returning their equivalents on the `to` form. Scores between tabulated
#' points are linearly interpolated; scores outside the range are held at the
#' nearest endpoint.
#'
#' @param eq An `equate_recipe` that has been through [run_equating()].
#' @param scores Numeric vector of raw scores on the `from` form.
#' @param from,to Optional form names selecting the pairing (see
#'   [conversion_table()]).
#' @param method Optional method name (e.g. `"Linear"`, `"Tucker"`,
#'   `"Equipercentile (No Smoothing)"`). If omitted and several are available,
#'   the first is used with a message.
#' @return A numeric vector of equated scores, the same length as `scores`.
#' @export
equated <- function(eq, scores, from = NULL, to = NULL, method = NULL) {
  conv <- eq@conversions
  if (is.null(conv) || !nrow(conv)) {
    cli::cli_abort("No conversions available. Run {.fn run_equating} first.")
  }
  sel <- .pick_pair(conv, from, to)
  sub <- conv[conv$from == sel$from & conv$to == sel$to, , drop = FALSE]

  available <- unique(sub$method)
  if (is.null(method)) {
    method <- available[1]
    if (length(available) > 1) {
      cli::cli_alert_info("Multiple methods; using {.val {method}}. Options: {.val {available}}.")
    }
  } else if (!method %in% available) {
    cli::cli_abort("Method {.val {method}} not found. Options: {.val {available}}.")
  }

  mm <- sub[sub$method == method, , drop = FALSE]
  mm <- mm[order(mm$x_score), ]
  stats::approxfun(mm$x_score, mm$equivalent_score, rule = 2)(scores)
}

#' Pull observed total scores for a pairing from the raw results
#' @keywords internal
.observed_scores <- function(eq, from, to) {
  grp <- eq@results[[paste0(from, ";", to)]]
  if (is.null(grp)) return(NULL)
  for (title in names(grp)) {
    g <- grp[[title]]
    if (!is.list(g) || inherits(g, "equate_failed")) next
    for (sm in names(g)) {
      r <- g[[sm]]
      if (is.list(r) && !is.null(r$observed_scores_x) && !is.null(r$observed_scores_y)) {
        return(list(x = as.numeric(r$observed_scores_x),
                    y = as.numeric(r$observed_scores_y)))
      }
    }
  }
  NULL
}

#' Moment comparison of observed and equated score distributions
#'
#' Rows: observed `from` scores, observed `to` scores, then each method's
#' equated `from` scores; columns: mean, sd, skewness, kurtosis.
#' @keywords internal
.conversion_moments <- function(eq, from, to) {
  ct  <- conversion_table(eq, from = from, to = to)
  obs <- .observed_scores(eq, from, to)
  mom <- function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE),
                       skewness = moments::skewness(x, na.rm = TRUE),
                       kurtosis = moments::kurtosis(x, na.rm = TRUE))

  rows <- list()
  if (!is.null(obs)) {
    rows[["Observed (from)"]] <- mom(obs$x)
    rows[["Observed (to)"]]   <- mom(obs$y)
  }
  for (m in setdiff(names(ct), "x_score")) {
    f <- stats::approxfun(ct$x_score, ct[[m]], rule = 2)
    vals <- if (!is.null(obs)) f(obs$x) else ct[[m]]
    rows[[m]] <- mom(vals)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
