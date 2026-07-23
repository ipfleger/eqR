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


# --- Numerical validation against independent reimplementations -----------------
# (self-contained; no external IRT package required)

.tcc <- function(a, b, cc, th) {
  vapply(th, function(t) sum(cc + (1 - cc) / (1 + exp(-a * (t - b)))), numeric(1))
}
.lw <- function(a, b, cc, th) {          # Lord-Wingersky: P(score = 0..n | theta)
  P <- 1
  for (i in seq_along(a)) {
    p <- cc[i] + (1 - cc[i]) / (1 + exp(-a[i] * (th - b[i])))
    P <- c(P * (1 - p), 0) + c(0, P * p)
  }
  P
}
.irt_fixture <- function() list(
  xit = paste0("X", 1:5), yit = paste0("Y", 1:5),
  ax = c(1.2, 0.8, 1.5, 1.0, 0.9), bx = c(-1.0, -0.3, 0.4, 1.1, 0.0), cx = c(.15, .2, .1, .18, .12),
  ay = c(1.0, 1.3, 0.7, 1.1, 1.4), by = c(-0.8, 0.1, 0.6, -0.2, 0.9), cy = c(.1, .15, .2, .12, .1),
  theta = seq(-4, 4, length.out = 100))
.irt_run <- function(f, params, type, ...) {
  init_equating() |>
    add_form(mk_form(300, f$xit, 0.5, seed = 1), "X", id_cols = "id", min_score = 0, max_score = 5, verbose = FALSE) |>
    add_form(mk_form(300, f$yit, 0.5, seed = 2), "Y", id_cols = "id", min_score = 0, max_score = 5, verbose = FALSE) |>
    add_plan(Y ~ X) |> add_design("RG") |>
    add_method("irt", irt_pars = params, type = type, theta = f$theta, ...) |>
    run_equating()
}

test_that("IRT true-score equating matches an independent uniroot TCC inversion", {
  f <- .irt_fixture()
  params <- data.frame(item = c(f$xit, f$yit), a = c(f$ax, f$ay), b = c(f$bx, f$by), c = c(f$cx, f$cy))
  pkg <- equated(.irt_run(f, params, "true_score"), 0:5, method = "IRT True Score")
  ind <- vapply(0:5, function(s) {
    lb <- sum(f$cx); ub <- length(f$ax)
    if (s <= lb || s >= ub) return(NA_real_)         # interior only (tails extrapolate)
    th <- uniroot(function(t) .tcc(f$ax, f$bx, f$cx, t) - s, c(-30, 30))$root
    .tcc(f$ay, f$by, f$cy, th)
  }, numeric(1))
  ok <- !is.na(ind)
  expect_lt(max(abs(pkg[ok] - ind[ok])), 1e-3)
  expect_true(all(diff(pkg) > 0))                    # strictly increasing
})

test_that("IRT observed-score equating matches an independent Lord-Wingersky computation", {
  f <- .irt_fixture()
  params <- data.frame(item = c(f$xit, f$yit), a = c(f$ax, f$ay), b = c(f$bx, f$by), c = c(f$cx, f$cy))
  pkg <- equated(.irt_run(f, params, "observed_score", theta_dist = "normal"), 0:5, method = "IRT Observed Score")
  w  <- dnorm(f$theta); w <- w / sum(w)
  dx <- rowSums(vapply(seq_along(f$theta), function(j) .lw(f$ax, f$bx, f$cx, f$theta[j]) * w[j], numeric(6)))
  dy <- rowSums(vapply(seq_along(f$theta), function(j) .lw(f$ay, f$by, f$cy, f$theta[j]) * w[j], numeric(6)))
  prd <- perc_rank(0:5, 0, 5, 1, cumsum(dx))
  ind <- EquiEquate(nsy = 6, miny = 0, incy = 1, crfdy = cumsum(dy), nsx = 6, prdx = prd)
  expect_lt(max(abs(pkg - ind)), 1e-3)
})

test_that("IRT equating of a form to itself is the identity", {
  f <- .irt_fixture()
  pid <- data.frame(item = c(f$xit, f$yit), a = c(f$ax, f$ax), b = c(f$bx, f$bx), c = c(f$cx, f$cx))
  ts <- equated(.irt_run(f, pid, "true_score"), 0:4, method = "IRT True Score")
  os <- equated(.irt_run(f, pid, "observed_score", theta_dist = "normal"), 0:5, method = "IRT Observed Score")
  expect_lt(max(abs(ts - 0:4)), 1e-3)   # interior; top score limited by the TCC ceiling
  expect_lt(max(abs(os - 0:5)), 1e-3)
})
