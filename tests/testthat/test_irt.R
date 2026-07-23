# IRT equating folds into the tidy result layer (no mirt needed: irt_coefs()
# accepts a plain data frame of item / a / b / c parameters).

mk_form <- function(n, items, p = 0.55, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::rbinom(n * length(items), 1, p), nrow = n,
              dimnames = list(NULL, items))
  cbind(data.frame(id = seq_len(n)), as.data.frame(m, check.names = FALSE))
}
mk_pars <- function(items, seed = 2) {
  set.seed(seed)
  data.frame(item = items,
             a = round(runif(length(items), 0.7, 1.8), 2),
             b = round(rnorm(length(items)), 2),
             c = 0)
}

test_that("RG IRT true- and observed-score equating fold into conversions()", {
  xitems <- paste0("xA", 1:10); yitems <- paste0("yB", 1:10)
  fx <- mk_form(500, xitems, 0.55, seed = 10)
  fy <- mk_form(520, yitems, 0.50, seed = 11)
  params <- mk_pars(c(xitems, yitems), seed = 12)  # concurrent calibration

  res <- init_equating() |>
    add_form(fx, name = "X", id_cols = "id", min_score = 0, max_score = 10, verbose = FALSE) |>
    add_form(fy, name = "Y", id_cols = "id", min_score = 0, max_score = 10, verbose = FALSE) |>
    add_plan(Y ~ X) |>
    add_design("RG") |>
    add_method("irt", irt_pars = params, type = "true_score") |>
    add_method("irt", irt_pars = params, type = "observed_score") |>
    run_equating()

  conv <- conversions(res)
  expect_true(all(c("IRT True Score", "IRT Observed Score") %in% conv$method))

  ct <- conversion_table(res)
  expect_equal(ct$x_score, 0:10)
  expect_true(all(diff(ct[["IRT True Score"]]) >= -1e-6))       # non-decreasing
  expect_true(all(diff(ct[["IRT Observed Score"]]) >= -1e-6))
  expect_length(equated(res, c(2, 5, 8), method = "IRT True Score"), 3)
})

test_that("IRT observed scores are full item totals (no dropped item)", {
  xitems <- paste0("xA", 1:10); yitems <- paste0("yB", 1:10)
  fx <- mk_form(300, xitems, 0.6, seed = 20)
  fy <- mk_form(300, yitems, 0.5, seed = 21)
  params <- mk_pars(c(xitems, yitems), seed = 22)

  res <- init_equating() |>
    add_form(fx, name = "X", id_cols = "id", min_score = 0, max_score = 10, verbose = FALSE) |>
    add_form(fy, name = "Y", id_cols = "id", min_score = 0, max_score = 10, verbose = FALSE) |>
    add_plan(Y ~ X) |>
    add_design("RG") |>
    add_method("irt", irt_pars = params, type = "true_score") |>
    run_equating()

  obs <- res@results[["X;Y"]][["R IRT true_score N"]][["IRT True Score"]]$observed_scores_x
  expect_length(obs, nrow(fx))
  expect_lte(max(obs), length(xitems))                 # would be 9 if an item were dropped
  expect_equal(obs, rowSums(fx[, xitems]))
})

test_that("common-item IRT runs the CNEG scale transformation and folds in", {
  common <- paste0("c", 1:5)
  xitems <- c(paste0("xU", 1:8), common)
  yitems <- c(paste0("yU", 1:8), common)
  fx <- mk_form(500, xitems, 0.55, seed = 30)
  fy <- mk_form(520, yitems, 0.50, seed = 31)
  px <- mk_pars(xitems, seed = 32); py <- mk_pars(yitems, seed = 33)

  res <- init_equating() |>
    add_form(fx, name = "X", id_cols = "id", min_score = 0, max_score = 13, verbose = FALSE) |>
    add_form(fy, name = "Y", id_cols = "id", min_score = 0, max_score = 13, verbose = FALSE) |>
    add_plan(Y ~ X) |>
    add_design("CG") |>
    add_method("irt", irt_pars = list(X = px, Y = py), type = "true_score") |>
    run_equating()

  conv <- conversions(res)
  expect_true("IRT True Score" %in% conv$method)
  ct <- conversion_table(res)
  expect_true(all(diff(ct[["IRT True Score"]]) >= -1e-6))
})
