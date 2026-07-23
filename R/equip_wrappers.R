# This file will contain the wrapper functions that will bridge my higher level syntax (e.g. run_equate) with the internal functions
# I am currently debugging


equipercentile <- function(forms, design, eq, title){
  if(design %in% c("S", "R")) {
    equipercentile_sgrg(eq = eq, forms = forms, title = title)
  } else if (design == "CG") {
    anchors <- get_anchors(eq)
    equipercentile_cg(eq = eq, forms = forms, anchors = anchors[[paste0(forms, collapse = ";")]], title = title)
  }
}

# --- Helper Functions for SG/RG Equipercentile Equating ---

#' @noRd
get_freq_dist <- function(scores, score_levels) {
  # This helper remains the same
  factor(scores, levels = score_levels) |> table() |> as.vector()
}

#' @noRd
equate_none_stat_fun <- function(data, i, score_ps, ...) {
  # data is now Long: cols are 'score', 'form'
  # i is the vector of indices selected by boot()

  sample_data <- data[i, ]

  # Extract X and Y scores based on the form column
  # We assume the forms are ordered as they were passed to the main function
  # (usually form X first, form Y second)

  # Identify the form labels from the score parameters or data
  # (Assuming the factor levels of 'form' match the order in score_ps)
  forms <- levels(sample_data$form)

  scores_x <- sample_data$score[sample_data$form == forms[1]]
  scores_y <- sample_data$score[sample_data$form == forms[2]]

  # Calculate frequency distributions
  freq_x <- get_freq_dist(scores_x, score_ps$range_x)
  freq_y <- get_freq_dist(scores_y, score_ps$range_y)

  # Proceed with equating as before
  crfd_x <- cumsum(freq_x / sum(freq_x))
  crfd_y <- cumsum(freq_y / sum(freq_y))

  prdx <- perc_rank(x = score_ps$range_x, min = score_ps$min_x, max = score_ps$max_x, inc = score_ps$inc_x, crfd = crfd_x)

  EquiEquate(nsy = length(score_ps$range_y),
             miny = score_ps$min_y,
             incy = score_ps$inc_y,
             crfdy = crfd_y,
             nsx = length(score_ps$range_x),
             prdx = prdx)
}

#' @noRd
equate_bb_stat_fun <- function(data, i, score_ps, method_opts) {
  sample_data <- data[i, ]

  # Extract scores by form
  forms <- levels(sample_data$form)
  scores_x <- sample_data$score[sample_data$form == forms[1]]
  scores_y <- sample_data$score[sample_data$form == forms[2]]

  freq_x <- get_freq_dist(scores_x, score_ps$range_x)
  freq_y <- get_freq_dist(scores_y, score_ps$range_y)

  # Smooth X
  smooth_x <- tryCatch(smooth_bb(n_persons = sum(freq_x), n_items = score_ps$k_x, freq = freq_x, rmoment = get_moments(scores=score_ps$range_x, freq=freq_x), nparm = method_opts$nparm, rel = score_ps$rel_x), error = function(e) NULL)
  crfd_x <- if (!is.null(smooth_x)) smooth_x$crfd else cumsum(freq_x / sum(freq_x))

  # Smooth Y
  smooth_y <- tryCatch(smooth_bb(n_persons = sum(freq_y), n_items = score_ps$k_y, freq = freq_y, rmoment = get_moments(scores=score_ps$range_y, freq=freq_y), nparm = method_opts$nparm, rel = score_ps$rel_y), error = function(e) NULL)
  crfd_y <- if (!is.null(smooth_y)) smooth_y$crfd else cumsum(freq_y / sum(freq_y))

  prdx <- perc_rank(x = score_ps$range_x, min = score_ps$min_x, max = score_ps$max_x, inc = score_ps$inc_x, crfd = crfd_x)
  EquiEquate(nsy = length(score_ps$range_y), miny = score_ps$min_y, incy = score_ps$inc_y, crfdy = crfd_y, nsx = length(score_ps$range_x), prdx = prdx)
}

#' @noRd
equate_loglinear_stat_fun <- function(data, i, score_ps, method_opts) {
  sample_data <- data[i, ]

  # Extract scores by form
  forms <- levels(sample_data$form)
  scores_x <- sample_data$score[sample_data$form == forms[1]]
  scores_y <- sample_data$score[sample_data$form == forms[2]]

  freq_x <- get_freq_dist(scores_x, score_ps$range_x)
  freq_y <- get_freq_dist(scores_y, score_ps$range_y)

  # Smooth X
  smooth_x <- tryCatch(smooth_ull(n = sum(freq_x), ns = length(freq_x), min = score_ps$min_x, inc = score_ps$inc_x, fd = freq_x, c = method_opts$degree), error = function(e) NULL)
  crfd_x <- if (!is.null(smooth_x) && smooth_x$converged) smooth_x$crfd else cumsum(freq_x / sum(freq_x))

  # Smooth Y
  smooth_y <- tryCatch(smooth_ull(n = sum(freq_y), ns = length(freq_y), min = score_ps$min_y, inc = score_ps$inc_y, fd = freq_y, c = method_opts$degree), error = function(e) NULL)
  crfd_y <- if (!is.null(smooth_y) && smooth_y$converged) smooth_y$crfd else cumsum(freq_y / sum(freq_y))

  prdx <- perc_rank(x = score_ps$range_x, min = score_ps$min_x, max = score_ps$max_x, inc = score_ps$inc_x, crfd = crfd_x)
  EquiEquate(nsy = length(score_ps$range_y), miny = score_ps$min_y, incy = score_ps$inc_y, crfdy = crfd_y, nsx = length(score_ps$range_x), prdx = prdx)
}

# Cubic-spline postsmoothing of equipercentile equating (symmetric).
# Follows Equating Recipes: smooth X->Y and Y->X separately with analytic
# equipercentile standard errors as weights, then report the average of the
# X->Y spline and the inverse of the Y->X spline. The spline is fit only over
# the score range whose percentile ranks fall within [prlow, prhigh].
#' @noRd
equate_spline_stat_fun <- function(data, i, score_ps, method_opts) {
  sp <- score_ps
  sd_   <- data[i, ]
  forms <- levels(sd_$form)
  sx <- sd_$score[sd_$form == forms[1]]
  sy <- sd_$score[sd_$form == forms[2]]
  nx <- length(sx); ny <- length(sy)

  fx <- get_freq_dist(sx, sp$range_x) / nx
  fy <- get_freq_dist(sy, sp$range_y) / ny
  crfdx <- cumsum(fx); crfdy <- cumsum(fy)
  prdx <- perc_rank(sp$range_x, sp$min_x, sp$max_x, sp$inc_x, crfdx)
  prdy <- perc_rank(sp$range_y, sp$min_y, sp$max_y, sp$inc_y, crfdy)

  win <- function(pr) {
    keep <- which(pr >= method_opts$prlow & pr <= method_opts$prhigh)
    if (length(keep) < 4L) keep <- seq_along(pr)   # fall back to the full range
    c(keep[1], keep[length(keep)])
  }

  # X -> Y direction
  eraw_xy <- EquiEquate(nsy = length(sp$range_y), miny = sp$min_y, incy = sp$inc_y,
                        crfdy = crfdy, nsx = length(sp$range_x), prdx = prdx)
  se_xy <- se_ep_equate(prdx, crfdy, nx, ny)
  wx <- win(prdx)
  dY <- post_smooth(sp$range_x, eraw_xy, se_xy, method_opts$s,
                    wx[1], wx[2], sp$max_y, sp$range_x)$vectY

  # Y -> X direction, then invert onto the X scale
  eraw_yx <- EquiEquate(nsy = length(sp$range_x), miny = sp$min_x, incy = sp$inc_x,
                        crfdy = crfdx, nsx = length(sp$range_y), prdx = prdy)
  se_yx <- se_ep_equate(prdy, crfdx, ny, nx)
  wy <- win(prdy)
  dX <- post_smooth(sp$range_y, eraw_yx, se_yx, method_opts$s,
                    wy[1], wy[2], sp$max_x, sp$range_y)$vectY
  dX_inv <- stats::approx(x = dX, y = sp$range_y, xout = sp$range_x, rule = 2)$y

  (dY + dX_inv) / 2
}

#' @noRd
equate_kernel_stat_fun <- function(data, i, score_ps, method_opts) {
  sample_data <- data[i, ]
  forms <- levels(sample_data$form)

  scores_x <- sample_data$score[sample_data$form == forms[1]]
  scores_y <- sample_data$score[sample_data$form == forms[2]]

  n_x <- length(scores_x)
  n_y <- length(scores_y)

  # Calculate relative frequencies per form
  rel_freq_x <- get_freq_dist(scores_x, score_ps$range_x) / n_x
  rel_freq_y <- get_freq_dist(scores_y, score_ps$range_y) / n_y

  # Calculate bandwidths (h) using the specific N for each form
  hx <- sd(scores_x) * (4 / (3 * n_x))^(1/5)
  hy <- sd(scores_y) * (4 / (3 * n_y))^(1/5)

  pr_x <- kernel_continu_cdf(score_ps$range_x, scores = score_ps$range_x, rel_freq = rel_freq_x, hx = hx)
  sapply(pr_x, function(p) kernel_inverse_cdf(p, scores = score_ps$range_y, rel_freq = rel_freq_y, hx = hy))
}

#' @noRd
equate_cll_sg_stat_fun <- function(data, i, score_ps = score_params, method_opts, cb = NULL) {
  sample_data <- data[i, ]
  forms <- levels(sample_data$form)

  scores_x <- sample_data$score[sample_data$form == forms[1]]
  scores_y <- sample_data$score[sample_data$form == forms[2]]

  # WARNING: Bivariate Log-Linear smoothing (CLL) requires paired data (Single Group design).
  # The "long" format with stratified bootstrapping breaks the pairing between X and Y.
  # This function will fail if N_x != N_y or if the data is not strictly paired.
  if (length(scores_x) != length(scores_y)) {
    stop("Bivariate smoothing requires equal sample sizes (paired data), but the input data is unpaired/unequal.")
  }

  # Attempt to create bivariate table assuming implicit pairing (row 1 x matches row 1 y)
  # Note: Stratified bootstrapping of the long format shuffles X and Y independently,
  # so the correlation structure here is likely destroyed.
  bdf <- table(factor(scores_x, levels = score_ps$range_x),
               factor(scores_y, levels = score_ps$range_y))

  bivar_smoothed <- smooth_bll(n = length(scores_x),
                               nsu = length(score_ps$range_x), minu = score_ps$min_x, incu = score_ps$inc_x,
                               nsv = length(score_ps$range_y), minv = score_ps$min_y, incv = score_ps$inc_y,
                               nct = as.vector(t(bdf)), scale = FALSE,
                               cu = method_opts$cu, cv = method_opts$cv, cuv = method_opts$cuv, cpm = method_opts$cpm)

  if (is.null(bivar_smoothed) || !bivar_smoothed$converged) stop("Bivariate log-linear smoothing failed to converge.")
  equate_cll(design = `if`(is.null(cb), "SG", "CB"), bivar = bivar_smoothed)
}

#' Master Dispatcher for SG/RG Equipercentile Statistic Functions
#' @noRd
equate_sgrg_statistic <- function(data, i, method_options, score_params, smooth_code, cb) {
  # This wrapper function selects the correct statistical engine based on the smoothing code.
  # It is designed to be passed as the `statistic` to boot::boot.
  switch(smooth_code,
         "N" = equate_none_stat_fun(data, i, score_ps = score_params, method_options),
         "B" = equate_bb_stat_fun(data, i, score_params, method_options),
         "L" = equate_loglinear_stat_fun(data, i, score_params, method_options),
         "S" = equate_spline_stat_fun(data, i, score_params, method_options), # error
         "K" = equate_kernel_stat_fun(data, i, score_params, method_options),
         "Z" = equate_cll_sg_stat_fun(data, i, score_params, method_options, cb = cb) # Fails to converge, doesn't have nct_b. I think we can fix this with our previous code.
  )
}


#' Robustly extract nested (50 and 95 percent) bootstrap CI bounds for one statistic
#'
#' Returns NA bounds for degenerate columns (all resamples equal) or when
#' `boot::boot.ci` fails (e.g. too few replications for "bca"), falling back
#' from the requested type to the percentile method. Prevents low-`boot_reps`
#' runs from aborting the whole equating.
#' @keywords internal
boot_ci_bounds <- function(bootobj, index, boot_type = "perc", conf = c(0.5, 0.95)) {
  na4 <- c(lower_bound_50 = NA_real_, upper_bound_50 = NA_real_,
           lower_bound_95 = NA_real_, upper_bound_95 = NA_real_)
  col <- bootobj$t[, index]
  if (length(unique(col[is.finite(col)])) < 2) return(na4)

  type <- boot_type[1]
  ci <- tryCatch(boot::boot.ci(bootobj, index = index, type = type, conf = conf),
                 error = function(e) NULL)
  if (is.null(ci) && type != "perc") {
    ci   <- tryCatch(boot::boot.ci(bootobj, index = index, type = "perc", conf = conf),
                     error = function(e) NULL)
    type <- "perc"
  }
  if (is.null(ci)) return(na4)

  comp <- switch(type, perc = "percent", norm = "normal",
                 basic = "basic", bca = "bca", "percent")
  m <- ci[[comp]]
  if (is.null(m) || !is.matrix(m) || nrow(m) < 2 || ncol(m) < 2) return(na4)
  nc <- ncol(m)
  c(lower_bound_50 = m[1, nc - 1], upper_bound_50 = m[1, nc],
    lower_bound_95 = m[2, nc - 1], upper_bound_95 = m[2, nc])
}

equipercentile_sgrg <- function(eq, forms, title) {

  # 1. Initial Setup
  method_options <- eq@methods[[title]]$options
  smooth_code <- eq@methods[[title]]$smooth
  method_name <- paste("Equipercentile", switch(smooth_code, N = "(No Smoothing)", B = "(Beta-Binomial)", L = "(Log-Linear)", S = "(Cubic Spline)", K = "(Kernel)", Z = "(CLL)"))

  boot_type <- method_options$boot_type
  boot_replications = method_options$boot_reps

  # 2. Data and Score Scale Preparation
  dat <- do.call(rbind, lapply(forms, function(frm) {
    data.frame(
      score = rowSums(eq@data[[frm]], na.rm = TRUE),
      form = frm,
      stringsAsFactors = FALSE # Ensure forms are treated consistently
    )
  }))
  # Ensure form is a factor for stratification
  dat$form <- factor(dat$form, levels = forms)

  score_params <- list(
    min_x = attr(eq@data[[forms[1]]], "min"), max_x = attr(eq@data[[forms[1]]], "max"), inc_x = attr(eq@data[[forms[1]]], "inc"), range_x = attr(eq@data[[forms[1]]], "range"), k_x = attr(eq@data[[forms[1]]], "k"), n_x = attr(eq@data[[forms[1]]], "n"), rel_x = attr(eq@data[[forms[1]]], "rel"),
    min_y = attr(eq@data[[forms[2]]], "min"), max_y = attr(eq@data[[forms[2]]], "max"), inc_y = attr(eq@data[[forms[2]]], "inc"), range_y = attr(eq@data[[forms[2]]], "range"), k_y = attr(eq@data[[forms[2]]], "k"), n_y = attr(eq@data[[forms[2]]], "n"), rel_y = attr(eq@data[[forms[2]]], "rel")
  )

  # 3. Bootstrapping and Result Formatting
  if (boot_replications <= 1) {
    equivalents <- equate_sgrg_statistic(data = dat,
                                         i = 1:nrow(dat),
                                         method_options = method_options,
                                         score_params = score_params, smooth_code = smooth_code,
                                         cb = `if`(is.null(attr(eq@design, "counter_balance_by")),
                                                   NULL,
                                                   eq@data_ids[[forms[1]]][[attr(eq@design, "counter_balance_by")]]))
    results <- list(list(
      parameters = NULL, x_score = score_params$range_x, equivalent_score = equivalents, bootstrapped_estimate = NA,
      nested_intervals = data.frame(se = NA, lower_bound_50 = NA, upper_bound_50 = NA, lower_bound_95 = NA, upper_bound_95 = NA),
      # MODIFIED: Filter observed scores from long format
      observed_scores_x = dat$score[dat$form == forms[1]],
      observed_scores_y = dat$score[dat$form == forms[2]]
    )) |> `names<-`(method_name)
  } else {
    equi_boot <- boot::boot(data = dat, R = boot_replications,
                            statistic = equate_sgrg_statistic,
                            # MODIFIED: Add stratification to resample within forms independently
                            strata = dat$form,
                            method_options = method_options,
                            score_params = score_params,
                            smooth_code = smooth_code,
                            cb = `if`(is.null(attr(eq@design, "counter_balance_by")),
                                      NULL,
                                      eq@data_ids[[forms[1]]][[attr(eq@design, "counter_balance_by")]]))

    cis <- do.call(rbind, lapply(seq_len(ncol(equi_boot$t)),
                                 function(i) boot_ci_bounds(equi_boot, i, boot_type)))

    bsm <- colMeans(equi_boot$t)
    bs_se <- apply(equi_boot$t, 2, sd)

    results <- list(list(
      parameters = NULL,
      x_score = attr(eq@data[[forms[1]]], "range"),
      equivalent_score = equi_boot$t0,
      bootstrapped_estimate = bsm,
      nested_intervals = data.frame(se = bs_se, cis),
      # MODIFIED: Filter observed scores from long format
      observed_scores_x = dat$score[dat$form == forms[1]],
      observed_scores_y = dat$score[dat$form == forms[2]]
    )) |> `names<-`(method_name)
  }

  return(results)
}

#' Prepare Data for Common-Item Equipercentile Equating
#'
#' @description
#' An internal helper function that transforms raw data from an equate_recipe
#' into the bivariate frequency distributions and score parameters needed for
#' CG equipercentile methods.
#'
#' @param eq The `equate_recipe` object.
#' @param forms A character vector of the two form names.
#' @param anchors A character vector of the anchor item names.
#'
#' @return A list containing bivariate frequency tables (`bdf_xv`, `bdf_yv`),
#'   score parameters (`score_params`), and sample sizes (`N1`, `N2`).
#' @noRd
prepare_cg_data <- function(eq, forms, anchors) {
  # Identify form names and data
  form_x_name <- forms[1]
  form_y_name <- forms[2]
  data_x <- eq@data[[form_x_name]]
  data_y <- eq@data[[form_y_name]]
  items_x <- eq@forms[[form_x_name]]
  items_y <- eq@forms[[form_y_name]]

  # Identify unique (non-anchor) items
  unique_items_x <- setdiff(items_x, anchors)
  unique_items_y <- setdiff(items_y, anchors)

  # Calculate total scores on unique items and anchor items for both populations
  scores_x <- rowSums(data_x[, unique_items_x, drop = FALSE], na.rm = TRUE)
  scores_v1 <- rowSums(data_x[, anchors, drop = FALSE], na.rm = TRUE)
  scores_y <- rowSums(data_y[, unique_items_y, drop = FALSE], na.rm = TRUE)
  scores_v2 <- rowSums(data_y[, anchors, drop = FALSE], na.rm = TRUE)

  # Determine the full, common score range for each variable
  min_x <- attr(eq@data[[form_x_name]], "min")
  max_x <- attr(eq@data[[form_x_name]], "max")
  x_range <- attr(eq@data[[form_x_name]], "range")

  min_y <- attr(eq@data[[form_y_name]], "min")
  max_y <- attr(eq@data[[form_y_name]], "max")
  y_range <-attr(eq@data[[form_y_name]], "range")

  min_v <- min(c(scores_v1, scores_v2, min_x, min_y))
  max_v <- sum(attr(eq@data[[form_x_name]], "points")[anchors])
  v_range <- seq(from = min_v, to = max_v, by = median(attr(eq@data[[form_x_name]], "points")[anchors]))

  # Create bivariate frequency distribution tables
  # Using factors ensures that all score points are represented, even with 0 frequency
  bdf_xv <- table(factor(scores_x, levels = x_range), factor(scores_v1, levels = v_range))
  bdf_yv <- table(factor(scores_y, levels = y_range), factor(scores_v2, levels = v_range))

  rel <- kr20(rbind(data_x[, anchors, drop = FALSE], data_y[, anchors, drop = FALSE]))
  # Assemble the list of score parameters
  score_params <- list(
    minx = min_x, maxx = max_x, incx = 1, nsx = length(x_range), rxv = rel,
    miny = min_y, maxy = max_y, incy = 1, nsy = length(y_range), ryv = rel,
    minv = min_v, maxv = max_v, incv = 1, nsv = length(v_range)
  )

  return(list(
    bdf_xv = bdf_xv,
    bdf_yv = bdf_yv,
    score_params = score_params,
    N1 = nrow(data_x),
    N2 = nrow(data_y)
  ))
}

#' Core Statistical Engine for CG Equipercentile Equating
#'
#' @description
#' This internal function is the core logic that is executed for each bootstrap
#' replication. It takes resampled frequency data, applies smoothing (if specified),
#' runs the appropriate equating methods, and returns the results.
#'
#' @param bdf_xv_boot,bdf_yv_boot Resampled bivariate frequency tables.
#' @param score_params A list of score scale parameters.
#' @param type The equating type code (e.g., 'E', 'C', 'A').
#' @param method_options A list of all method options.
#'
#' @return A named list of equated score vectors.
#' @noRd
equate_cg_statistic <- function(bdf_xv_boot, bdf_yv_boot, score_params, type, method_options) {

  # Extract score parameters for convenience
  sp <- score_params

  # Calculate marginal distributions from the bootstrapped bivariate tables
  fx_v1 <- rowSums(bdf_xv_boot)
  gv_v1 <- colSums(bdf_xv_boot)
  gy_v2 <- rowSums(bdf_yv_boot)
  gv_v2 <- colSums(bdf_yv_boot)

  # --- Apply Smoothing (Future Enhancement) ---
  # For now, we only implement the "No Smoothing" case.
  # A switch statement would go here to handle different smooth_codes.
  crfd_x1 <- cumsum(fx_v1 / sum(fx_v1))
  crfd_y2 <- cumsum(gy_v2 / sum(gy_v2))
  crfd_v1 <- cumsum(gv_v1 / sum(gv_v1))
  crfd_v2 <- cumsum(gv_v2 / sum(gv_v2))

  # --- Execute Equating Method(s) ---
  results <- list()

  # Helper to run Frequency Estimation / MFE
  run_fe <- function(is_mfe = FALSE) {
    rv1 <- ifelse(is_mfe, method_options$rel %||% sp$rxv, 0)
    rv2 <- ifelse(is_mfe, method_options$rel %||% sp$ryv, 0)
    fe_mfe_equate(w1 = method_options$w1, internal = method_options$internal_anchors,
                  bxvin = bdf_xv_boot, byvin = bdf_yv_boot,
                  rv1 = rv1, rv2 = rv2,
                  minx = sp$minx, maxx = sp$maxx, incx = sp$incx,
                  miny = sp$miny, maxy = sp$maxy, incy = sp$incy)
  }

  # Helper to run Chained Equipercentile
  run_chained <- function() {
    prdx1 <- perc_rank(x = sp$minx:sp$maxx, min = sp$minx, max = sp$maxx, inc = sp$incx, crfd = crfd_x1)
    chained_equate(nsx = sp$nsx, prx1 = prdx1,
                   minv = sp$minv, maxv = sp$maxv, incv = sp$incv, nsv = sp$nsv, Gv1 = crfd_v1,
                   miny = sp$miny, incy = sp$incy, nsy = sp$nsy, Gy2 = crfd_y2, Gv2 = crfd_v2)
  }

  # Use a switch-like structure for the equating type
  if (type %in% c("E", "G", "H", "A")) {
    fe_results <- run_fe()
    results$FrequencyEstimation <- fe_results$eraw
    results$BraunHolland_FE <- fe_results$eraw_bh
  }
  if (type %in% c("F", "G", "A")) {
    mfe_results <- run_fe(is_mfe = TRUE)
    results$ModifiedFrequencyEstimation <- mfe_results$eraw
    results$BraunHolland_MFE <- mfe_results$eraw_bh
  }
  if (type %in% c("C", "H", "A")) {
    results$Chained <- run_chained()
  }

  return(results)
}


equipercentile_cg <- function(eq, forms, anchors, title) {#, boot_type = "perc", boot_replications = 1000

  # 1. Initial Setup
  method_spec <- eq@methods[[title]]
  method_options <- eq@methods[[title]]$options
  type <- method_spec$type

  boot_type <- method_options$boot_type
  boot_replications = method_options$boot_reps

  # 2. Prepare Data
  prep_data <- prepare_cg_data(eq, forms, anchors)
  score_params <- prep_data$score_params

  # Get original scores for the summary object
  dat <- data.frame(do.call(cbind, lapply(forms |> `names<-`(forms), \(frm){
    rowSums(eq@data[[frm]][eq@forms[[frm]]], na.rm = TRUE)
  })))
  names(dat) <- c("x", "y")

  # 3. Pre-run to get point estimates and structure
  point_estimates_list <- equate_cg_statistic(bdf_xv_boot = prep_data$bdf_xv,
                                              bdf_yv_boot = prep_data$bdf_yv,
                                              score_params = score_params,
                                              type = type,
                                              method_options = method_options)
  method_names <- names(point_estimates_list)
  n_scores <- length(point_estimates_list[[1]])

  # 4. Bootstrap or format results
  if (boot_replications <= 1) {
    results <- lapply(method_names, function(name) {
      list(
        parameters = NULL,
        x_score = score_params$minx:score_params$maxx,
        equivalent_score = point_estimates_list[[name]],
        bootstrapped_estimate = NA,
        nested_intervals = data.frame(se = NA, lower_bound_50 = NA, upper_bound_50 = NA, lower_bound_95 = NA, upper_bound_95 = NA),
        observed_scores_x = dat$x,
        observed_scores_y = dat$y
      )
    }) |> `names<-`(method_names)
  } else {
    # Define the statistic function for boot::boot
    boot_stat_fun <- function(data, i) {
      # Resample the bivariate frequency distributions
      # The 'data' argument is a placeholder; we use the indices 'i' to resample our original tables
      indices1 <- i[1:prep_data$N1]
      indices2 <- i[(prep_data$N1 + 1):(prep_data$N1 + prep_data$N2)]

      # This is a conceptual resampling; a real implementation would resample the raw scores
      # and rebuild the tables. For simplicity with multinomial resampling of tables:
      bdf_xv_boot <- bootstrap_bivariate(prep_data$N1, prep_data$bdf_xv)
      bdf_yv_boot <- bootstrap_bivariate(prep_data$N2, prep_data$bdf_yv)

      # Run the equating engine
      boot_results_list <- equate_cg_statistic(bdf_xv_boot, bdf_yv_boot, score_params, type, method_options)

      # Unlist to return a single vector
      unlist(boot_results_list, use.names = FALSE)
    }

    # We need a dummy data object for boot to get the sample size right
    dummy_data <- 1:(prep_data$N1 + prep_data$N2)
    equi_boot <- boot::boot(data = dummy_data, statistic = boot_stat_fun, R = boot_replications,
                            strata = c(rep(1, prep_data$N1), rep(2, prep_data$N2)))

    # Process bootstrap results
    bsm <- colMeans(equi_boot$t)
    bs_se <- apply(equi_boot$t, 2, sd)

    results <- lapply(seq_along(method_names), function(m_idx) {
      name <- method_names[m_idx]
      # Define column indices for this method in the bootstrap output matrix
      start_col <- (m_idx - 1) * n_scores + 1
      end_col <- m_idx * n_scores
      cols <- start_col:end_col

      # Calculate CIs for this method's columns
      cis <- do.call(rbind, lapply(cols, function(i) boot_ci_bounds(equi_boot, i, boot_type)))

      list(
        parameters = NULL,
        x_score = score_params$minx:score_params$maxx,
        equivalent_score = equi_boot$t0[cols],
        bootstrapped_estimate = bsm[cols],
        nested_intervals = data.frame(se = bs_se[cols], cis),
        observed_scores_x = dat$x,
        observed_scores_y = dat$y
      )
    }) |> `names<-`(method_names)
  }

  return(results)
}
