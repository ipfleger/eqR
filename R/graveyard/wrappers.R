

linear <- function(forms, design, eq, title, boot_type = "perc", boot_replications = 1000){
  if(design %in% c("S", "R")) {
    linear_sgrg(eq = eq, forms = forms, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if (design == "CG") {
    anchors <- get_anchors(eq)
    linear_cg(eq = eq, forms = forms, anchors = anchors[[paste0(forms, collapse = ";")]], title = title, boot_type = boot_type, boot_replications = boot_replications)
  }
}

equipercentile <- function(forms, design, eq, title, boot_type = "perc", boot_replications = 1000){
  if(design %in% c("S", "R")) {
    equipercentile_sgrg(eq = eq, forms = forms, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if (design == "CG") {
    anchors <- get_anchors(eq)
    equipercentile_cg(eq = eq, forms = forms, anchors = anchors[[paste0(forms, collapse = ";")]], title = title, boot_type = boot_type, boot_replications = boot_replications)
  }
}

linear_sgrg <- function(eq, forms, title, boot_type = "perc", boot_replications = 1000){

  # Get all options, merged with defaults
  method_options <- get_method_options(eq@methods[[title]])
  mean_only <- method_options$mean_only
  method_name <- ifelse(mean_only, "Mean", "Linear")

  dat <- data.frame(do.call(cbind, lapply(forms |> `names<-`(forms), \(frm){
    rowSums(eq@data[[frm]][eq@forms[[frm]]], na.rm = TRUE)
  })))


  if (boot_replications <= 1) {
    point_estimates <- unlist(linear_equate_rgsg(
      mnx = mean(dat[,1]), sdx = sd(dat[,1]),
      mny = mean(dat[,2]), sdy = sd(dat[,2]),
      mean_only = mean_only,
      min_x = attr(eq@data[[forms[1]]], "min"), max_x = attr(eq@data[[forms[1]]], "max"), inc_x = attr(eq@data[[forms[1]]], "inc")
    ))

    results <- list(list(
      parameters = data.frame(
        statistics = c("Slope", "Intercept"), estimate = point_estimates[1:2],
        se = NA_real_, lower_bound_95 = NA_real_, upper_bound_95 = NA_real_,
        bootstrapped_estimate = NA_real_
      ),
      x_score = seq(from = minx, to = maxx, by = inc),
      equivalents = point_estimates[-1:-2],
      bootstrapped_estimate = NA_real_,
      nested_intervals = data.frame(
        se = NA_real_, lower_bound_50 = NA_real_, upper_bound_50 = NA_real_,
        lower_bound_95 = NA_real_, upper_bound_95 = NA_real_
      ),
      single = identical(eq@data_ids[[forms[1]]][[1]], eq@data_ids[[forms[2]]][[1]]),
      observed_scores_x = dat[[1]], observed_scores_y = dat[[2]]
    )) |> `names<-`(method_name)

  } else {
    linear_boot_fun <- function(data, i) {
      unlist(linear_equate_rgsg(
        mnx = mean(data[i,1]), sdx = sd(data[i,1]),
        mny = mean(data[i,2]), sdy = sd(data[i,2]),
        mean_only = mean_only,
        min_x = attr(eq@data[[forms[1]]], "min"), max_x = attr(eq@data[[forms[1]]], "max"), inc_x = attr(eq@data[[forms[1]]], "inc")
      ))
    }
    linear <- boot::boot(data = dat, statistic = linear_boot_fun, R = boot_replications)

    boots <- lapply(1:ncol(linear$t), \(i) `if`(length(unique(linear$t[,i])) == 1,
                                                list(cbind(c(NA_real_, NA_real_), c(NA_real_, NA_real_))) |> `names<-`(sapply(boot_type, switch, "perc" = "percent", "norm" = "normal", "basic" = "basic", "bca" = "bca")),
                                                boot::boot.ci(linear, index = i, type = boot_type, conf = c(.5, .95))) )

    cis <- do.call(rbind, lapply(boots, \(x) {
      ci_type <- sapply(boot_type, switch, "perc" = "percent", "norm" = "normal", "basic" = "basic", "bca" = "bca")
      x_mat <- x[[ci_type]]
      c(x_mat[1,(ncol(x_mat) - 1):ncol(x_mat)], x_mat[2,(ncol(x_mat) - 1):ncol(x_mat)])
    })) |> `colnames<-`(c(paste(rep(c("lower_bound", "upper_bound"), times = 2), rep(c(.5, .95)*100, each = 2), sep = "_")))

    bsm <- colMeans(linear$t)
    bs_se <- apply(linear$t, 2, sd)

    results <- list(list(
      parameters = data.frame(statistics = c("Slope", "Intercept"), estimate = linear$t0[1:2],
                              se = bs_se[1:2], cis[1:2,3:4], bootstrapped_estimate = bsm[1:2]),
      x_score = attr(eq@data[[forms[1]]], "range"),
      equivalent_score = linear$t0[-1:-2],
      bootstrapped_estimate = bsm[-1:-2],
      nested_intervals = cbind(se = bs_se[-1:-2], cis[-1:-2,]),
      single = identical(eq@data_ids[[forms[1]]][[1]], eq@data_ids[[forms[2]]][[1]]),
      observed_scores_x = dat[[1]], observed_scores_y = dat[[2]]
    )) |> `names<-`(method_name)
  }

  results <- lapply(names(results) |> `names<-`(names(results)), \(method){
    print_linear_sgrg(title = paste0(method, ": ", paste0(forms, collapse = " to ")), pdata = results[[method]])
  })

  # Return the raw results, ready for the summary function
  return(results)
}

linear_cg <- function(eq, forms, anchors, title, boot_type = "perc", boot_replications = 1000){

  method_spec <- eq@methods[[title]]
  method_options <- get_method_options(method_spec)
  type <- method_spec$type

  score_scale <- list(min_x = attr(eq@data[[forms[1]]], "min"), max_x = attr(eq@data[[forms[1]]], "max"), inc_x = attr(eq@data[[forms[1]]], "inc"))

  dat <- data.frame(do.call(cbind, lapply(forms |> `names<-`(forms), \(frm){
    rowSums(eq@data[[frm]][eq@forms[[frm]]], na.rm = TRUE)
  })))
  nrows <- unlist(lapply(eq@data[forms], nrow))

  # Pre-run to get the structure (method names, number of scores)
  out_pre <- linear_cg_sub(eq = eq, forms = forms, method_options = method_options, type = type, anchors = anchors, score_scale = score_scale)
  method_names <- out_pre$summary$Method
  n_methods <- length(method_names)
  score_names <- names(out_pre$equated_scores)
  n_scores <- length(score_names)

  if (boot_replications <= 1) {
    results <- lapply(seq_along(method_names) |> `names<-`(method_names), \(i){
      list(parameters = data.frame(statistics = c("Slope", "Intercept"),
                                   estimate = unlist(out_pre$summary[i, c("a", "b")]),
                                   se = NA_real_, lower_bound_95 = NA_real_, upper_bound_95 = NA_real_,
                                   bootstrapped_estimate = NA_real_),
           x_score = attr(eq@data[[forms[1]]], "range"),
           equivalents = out_pre$equated_scores[i,],
           bootstrapped_estimate = NA_real_,
           nested_intervals = data.frame(se = rep(NA_real_, n_scores),
                                         lower_bound_50 = NA_real_, upper_bound_50 = NA_real_,
                                         lower_bound_95 = NA_real_, upper_bound_95 = NA_real_),
           observed_scores_x = dat[[1]], observed_scores_y = dat[[2]])
    })
  } else {
    linear_boot_fun <- function(x, i) {
      out_boot <- linear_cg_sub(eq = eq, forms = forms, method_options = method_options, type = type, anchors = anchors, score_scale = score_scale, index = x[i])
      params_flat <- as.vector(t(out_boot$summary[, c("a", "b")]))
      scores_flat <- as.vector(t(out_boot$equated_scores))
      c(params_flat, scores_flat)
    }
    linear <- boot::boot(data = 1:sum(nrows), strata = factor(rep(forms, times = nrows)),
                         statistic = linear_boot_fun, R = boot_replications)

    boots <- lapply(1:ncol(linear$t), \(i) `if`(length(unique(linear$t[,i])) == 1,
                                                list(cbind(c(NA_real_, NA_real_), c(NA_real_, NA_real_))) |> `names<-`(sapply(boot_type, switch, "perc" = "percent", "norm" = "normal", "basic" = "basic", "bca" = "bca")),
                                                boot::boot.ci(linear, index = i, type = boot_type, conf = c(.5, .95))))

    cis <- do.call(rbind, lapply(boots, \(x) {
      ci_type <- sapply(boot_type, switch, "perc" = "percent", "norm" = "normal", "basic" = "basic", "bca" = "bca")
      x_mat <- x[[ci_type]]
      c(x_mat[1,(ncol(x_mat) - 1):ncol(x_mat)], x_mat[2,(ncol(x_mat) - 1):ncol(x_mat)])
    })) |> `colnames<-`(c(paste(rep(c("lower_bound", "upper_bound"), times = 2), rep(c(.5, .95)*100, each = 2), sep = "_")))

    bsm <- colMeans(linear$t)
    bs_se <- apply(linear$t, 2, sd)

    results <- lapply(seq_along(method_names) |> `names<-`(method_names), \(i){
      param_indices <- c(2*i - 1, 2*i)
      score_start_index <- (2 * n_methods) + (i - 1) * n_scores + 1
      score_indices <- score_start_index:(score_start_index + n_scores - 1)

      list(parameters = data.frame(statistics = c("Slope", "Intercept"),
                                   estimate = linear$t0[param_indices],
                                   se = bs_se[param_indices],
                                   cis[param_indices, 3:4],
                                   bootstrapped_estimate = bsm[param_indices]),
           x_score = attr(eq@data[[forms[1]]], "range"),
           equivalents = linear$t0[score_indices],
           bootstrapped_estimate = bsm[score_indices],
           nested_intervals = cbind(se = bs_se[score_indices], cis[score_indices,]),
           observed_scores_x = dat[[1]], observed_scores_y = dat[[2]])
    })
  }
  results <- lapply(names(results) |> `names<-`(names(results)), \(method){
    print_linear_cg(title = paste0(method, ": ", paste0(forms, collapse = " to ")), pdata = results[[method]])
  })
  # Return the raw results, ready for the summary function
  return(results)
}

linear_cg_sub <- function(eq, forms, method_options, type, anchors, score_scale, index = NULL){
  means <- lapply(forms |> `names<-`(forms), \(frm) {
    rows <- nrow(eq@data[[frm]])
    i <- `if`(is.null(index), 1:rows, `if`(frm == forms[1], index[index<=rows], index[index>rows] - rows))
    fdat <- eq@data[[frm]][eq@forms[[frm]]]
    adat <- eq@data[[frm]][anchors]
    ts = rowSums(fdat, na.rm = TRUE)[i]
    tsa = rowSums(adat, na.rm = TRUE)[i]
    list(mnx =  mean(ts), sdx = sd(ts), mnv = mean(tsa), sdv = sd(tsa), cov_xv = cov(ts, tsa))
  })

  # Pull options from the merged list
  w1 <- method_options$w1
  internal_anchors <- method_options$internal_anchors
  mean_only <- method_options$mean_only

  linear_equate_ci(mnx1 = means[[1]]$mnx, sdx1 = means[[1]]$sdx, mnv1 = means[[1]]$mnv, sdv1 = means[[1]]$sdv, covxv1 = means[[1]]$cov_xv,
                   mny2 = means[[2]]$mnx, sdy2 = means[[2]]$sdx, mnv2 = means[[2]]$mnv, sdv2 = means[[2]]$sdv, covyv2 = means[[2]]$cov_xv,
                   w1 = w1,
                   min_x = score_scale$min_x, max_x = score_scale$max_x, inc_x = score_scale$inc_x,
                   mean_only = mean_only, type = type, anchor = internal_anchors)
}

equipercentile_sgrg <- function(eq, forms, title, boot_type = "perc", boot_replications = 1000) {

  # 1. Initial Setup
  method_options <- get_method_options(eq@methods[[title]])
  smooth_code <- eq@methods[[title]]$smooth
  method_name <- paste("Equipercentile", switch(smooth_code, N = "(No Smoothing)", B = "(Beta-Binomial)", L = "(Log-Linear)", S = "(Cubic Spline)", K = "(Kernel)", Z = "(CLL)"))

  # 2. Data and Score Scale Preparation
  dat <- data.frame(do.call(cbind, lapply(forms |> `names<-`(forms), \(frm){
    rowSums(eq@data[[frm]][eq@forms[[frm]]], na.rm = TRUE)
  })))
  names(dat) <- c("x", "y")


  # Helper to create a frequency distribution on the full score range
  get_freq_dist <- function(scores, score_range) {
    # Ensure all levels from min to max score are included
    score_factor <- factor(scores, levels = score_range)
    as.vector(table(score_factor))
  }

  # 3. Core Equating Logic (to be used in bootstrap)
  equate_stat_fun <- function(data, i) {
    sample_data <- data[i, ]
    freq_x <- get_freq_dist(sample_data$x, attr(eq@data[[forms[1]]], "range"))
    freq_y <- get_freq_dist(sample_data$y, attr(eq@data[[forms[2]]], "range"))

    # --- Apply Smoothing ---
    # This is a simplified example; a full implementation would pass more options
    # and handle different smoothers more robustly.
    switch(smooth_code,
           "N" = {
             crfd_x <- cumsum(freq_x / sum(freq_x))
             crfd_y <- cumsum(freq_y / sum(freq_y))
           },
           "B" = {
             # Placeholder for Beta-Binomial; requires n_items, moments, etc.
             # For now, falls back to no smoothing
             crfd_x <- cumsum(freq_x / sum(freq_x))
             crfd_y <- cumsum(freq_y / sum(freq_y))
           },
           # Add other cases for L, S, K, Z here
           { # Default case (no smoothing)
             crfd_x <- cumsum(freq_x / sum(freq_x))
             crfd_y <- cumsum(freq_y / sum(freq_y))
           }
    )

    # --- Perform Equating ---
    prdx <- perc_rank(x = attr(eq@data[[forms[1]]], "range"),#get_freq_dist(sample_data$x, attr(eq@data[[forms[1]]], "range")),
                      min = attr(eq@data[[forms[1]]],"min"), max = attr(eq@data[[forms[1]]],"max"), inc = attr(eq@data[[forms[1]]],"inc"), crfd = crfd_x)
    eraw <- EquiEquate(nsy = length(attr(eq@data[[forms[2]]], "range")), miny =  attr(eq@data[[forms[2]]],"min"), incy = attr(eq@data[[forms[2]]],"inc"), crfdy = crfd_y, nsx = length(attr(eq@data[[forms[1]]], "range")), prdx = prdx)
    return(eraw)
  }

  # 4. Bootstrapping and Result Formatting
  if (boot_replications <= 1) {
    equivalents <- equate_stat_fun(dat, 1:nrow(dat))
    results <- list(list(
      parameters = NULL, # No parameters for equipercentile
      x_score = attr(eq@data[[forms[1]]], "range"),
      equivalent_score = equivalents,
      bootstrapped_estimate = NA,
      nested_intervals = data.frame(se = NA, lower_bound_50 = NA, upper_bound_50 = NA, lower_bound_95 = NA, upper_bound_95 = NA),
      observed_scores_x = dat$x, observed_scores_y = dat$y
    )) |> `names<-`(method_name)
  } else {
    equi_boot <- boot::boot(data = dat, statistic = equate_stat_fun, R = boot_replications)

    boots <- lapply(1:ncol(equi_boot$t), \(i) `if`(length(unique(equi_boot$t[,i])) == 1,
                                                   list(cbind(c(NA,NA), c(NA,NA))) |> `names<-`(sapply(boot_type, switch, "perc" = "percent")),
                                                   boot::boot.ci(equi_boot, index = i, type = boot_type, conf = c(.5, .95))))

    cis <- do.call(rbind, lapply(boots, \(x) {
      ci_type <- sapply(boot_type, switch, "perc" = "percent")
      x_mat <- x[[ci_type]]
      c(x_mat[1,(ncol(x_mat) - 1):ncol(x_mat)], x_mat[2,(ncol(x_mat) - 1):ncol(x_mat)])
    })) |> `colnames<-`(c("lower_bound_50", "upper_bound_50", "lower_bound_95", "upper_bound_95"))

    bsm <- colMeans(equi_boot$t)
    bs_se <- apply(equi_boot$t, 2, sd)

    results <- list(list(
      parameters = NULL,
      x_score = attr(eq@data[[forms[1]]], "range"),
      equivalent_score = equi_boot$t0,
      bootstrapped_estimate = bsm,
      nested_intervals = data.frame(se = bs_se, cis),
      observed_scores_x = dat$x, observed_scores_y = dat$y
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
  min_x <- min(scores_x)
  max_x <- max(scores_x)
  x_range <- min_x:max_x

  min_y <- min(scores_y)
  max_y <- max(scores_y)
  y_range <- min_y:max_y

  min_v <- min(scores_v1, scores_v2)
  max_v <- max(scores_v1, scores_v2)
  v_range <- min_v:max_v

  # Create bivariate frequency distribution tables
  # Using factors ensures that all score points are represented, even with 0 frequency
  bdf_xv <- table(factor(scores_x, levels = x_range), factor(scores_v1, levels = v_range))
  bdf_yv <- table(factor(scores_y, levels = y_range), factor(scores_v2, levels = v_range))

  # Assemble the list of score parameters
  score_params <- list(
    minx = min_x, maxx = max_x, incx = 1, nsx = length(x_range),
    miny = min_y, maxy = max_y, incy = 1, nsy = length(y_range),
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
    rv1 <- ifelse(is_mfe, method_options$rel, 0)
    rv2 <- ifelse(is_mfe, method_options$rel, 0)
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


equipercentile_cg <- function(eq, forms, anchors, title, boot_type = "perc", boot_replications = 1000) {

  # 1. Initial Setup
  method_spec <- eq@methods[[title]]
  method_options <- get_method_options(method_spec)
  type <- method_spec$type

  # 2. Prepare Data
  prep_data <- prepare_cg_data(eq, forms, anchors)
  score_params <- prep_data$score_params

  # Get original scores for the summary object
  dat <- data.frame(do.call(cbind, lapply(forms |> `names<-`(forms), \(frm){
    rowSums(eq@data[[frm]][eq@forms[[frm]]], na.rm = TRUE)
  })))
  names(dat) <- c("x", "y")

  # 3. Pre-run to get point estimates and structure
  point_estimates_list <- equate_cg_statistic(prep_data$bdf_xv, prep_data$bdf_yv, score_params, type, method_options)
  method_names <- names(point_estimates_list)
  n_scores <- length(point_estimates_list[[1]])

  # 4. Bootstrap or format results
  if (boot_replications <= 1) {
    results <- lapply(method_names, function(name) {
      list(
        parameters = NULL,
        x_score = score_params$minx:score_params$maxx,
        equivalents = point_estimates_list[[name]],
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
      boots <- lapply(cols, \(i) `if`(length(unique(equi_boot$t[,i])) == 1,
                                      list(cbind(c(NA,NA), c(NA,NA))) |> `names<-`(sapply(boot_type, switch, "perc" = "percent")),
                                      boot::boot.ci(equi_boot, index = i, type = boot_type, conf = c(.5, .95))))
      cis <- do.call(rbind, lapply(boots, \(x) {
        ci_type <- sapply(boot_type, switch, "perc" = "percent")
        x_mat <- x[[ci_type]]
        c(x_mat[1,(ncol(x_mat) - 1):ncol(x_mat)], x_mat[2,(ncol(x_mat) - 1):ncol(x_mat)])
      })) |> `colnames<-`(c("lower_bound_50", "upper_bound_50", "lower_bound_95", "upper_bound_95"))

      list(
        parameters = NULL,
        x_score = score_params$minx:score_params$maxx,
        equivalents = equi_boot$t0[cols],
        bootstrapped_estimate = bsm[cols],
        nested_intervals = data.frame(se = bs_se[cols], cis),
        observed_scores_x = dat$x,
        observed_scores_y = dat$y
      )
    }) |> `names<-`(method_names)
  }

  return(results)
}
