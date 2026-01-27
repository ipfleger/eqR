

linear <- function(forms, design, eq, title){
  if(design %in% c("S", "R")) {
    linear_sgrg(eq = eq, forms = forms, title = title)
  } else if (design == "CG") {
    anchors <- get_anchors(eq)
    linear_cg(eq = eq, forms = forms, anchors = anchors[[paste0(forms, collapse = ";")]], title = title)
  }
}

linear_sgrg <- function(eq, forms, title){

  # Get all options, merged with defaults
  method_options <- eq@methods[[title]]$options
  mean_only <- method_options$mean_only
  method_name <- ifelse(mean_only, "Mean", "Linear")

  boot_type <- method_options$boot_type
  boot_replications = method_options$boot_reps

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
      equivalent_score = point_estimates[-1:-2],
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

  # results <- lapply(names(results) |> `names<-`(names(results)), \(method){
  #   print_linear_sgrg(title = paste0(method, ": ", paste0(forms, collapse = " to ")), pdata = results[[method]])
  # })

  # Return the raw results, ready for the summary function
  return(results) # Someday we might allow multiple, but probably not I guess...
}

#' Perform Linear Equating for Random Groups and Single Group Designs
#'
#' @description
#' Performs linear or mean equating for the Random Groups (RG) and Single Group
#' (SG) designs. This is an R translation of the `RGandSG_LinEq` function.
#'
#' @details
#' The function calculates the slope (`a`) and intercept (`b`) of the linear
#' transformation `y = a*x + b`. For mean equating, the slope is fixed at 1.
#' For linear equating, the slope is the ratio of the standard deviations.
#' It then applies this transformation to a sequence of raw scores.
#'
#' @param mnx A numeric value, the mean of the new form (X).
#' @param sdx A numeric value, the standard deviation of the new form (X).
#' @param mny A numeric value, the mean of the old form (Y).
#' @param sdy A numeric value, the standard deviation of the old form (Y).
#' @param mean_only A logical. If TRUE, mean equating is performed (slope is fixed to 1).
#'   If FALSE (the default), linear equating is performed.
#' @param min_x A numeric value, the minimum score on the X scale.
#' @param max_x A numeric value, the maximum score on the X scale.
#' @param inc_x A numeric value, the increment between scores on the X scale.
#'
#' @return A list containing the slope (`a`), intercept (`b`), and a numeric
#'   vector of the equated raw scores (`equated_scores`).
#' @author R. L. Brennan (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' \dontrun{
#' # Example data for linear equating
#' mean_x <- 25; sd_x <- 5
#' mean_y <- 50; sd_y <- 10
#' min_score <- 10; max_score <- 40; increment <- 1
#'
#' result <- linear_equate_rgsg(mean_x, sd_x, mean_y, sd_y, mean_only = FALSE,
#'                              min_score, max_score, increment)
#' print(result$a) # Slope should be 10/5 = 2
#' print(result$b) # Intercept should be 50 - 2*25 = 0
#' print(head(result$equated_scores))
#' }
linear_equate_rgsg <- function(mnx, sdx, mny, sdy, mean_only = FALSE, min_x, max_x, inc_x) {

  # Determine the slope 'a' based on the method
  if (mean_only) {
    a <- 1.0
  } else {
    if (sdx == 0) {
      stop("Standard deviation of X (sdx) cannot be zero for linear equating. use option mean_only = TRUE for mean equating.")
    }
    a <- sdy / sdx
  }

  # Calculate the intercept 'b'
  b <- mny - a * mnx

  # Generate the sequence of raw scores for form X
  raw_scores_x <- seq(from = min_x, to = max_x, by = inc_x)

  # Apply the linear transformation to get equated scores
  equated_scores <- b + a * raw_scores_x

  # Return the results in a list
  return(list(a = a, b = b, equated_scores = equated_scores))
}


linear_cg <- function(eq, forms, anchors, title){

  method_spec <- eq@methods[[title]]
  method_options <- eq@methods[[title]]$options
  type <- method_spec$type

  boot_type <- method_options$boot_type
  boot_replications = method_options$boot_reps


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
           equivalent_score = out_pre$equated_scores[i,],
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
           equivalent_score = linear$t0[score_indices],
           bootstrapped_estimate = bsm[score_indices],
           nested_intervals = cbind(se = bs_se[score_indices], cis[score_indices,]),
           observed_scores_x = dat[[1]], observed_scores_y = dat[[2]])
    })
  }
  # results <- lapply(names(results) |> `names<-`(names(results)), \(method){
  #   print_linear_cg(title = paste0(method, ": ", paste0(forms, collapse = " to ")), pdata = results[[method]])
  # })
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

