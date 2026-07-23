# Tests for the tidy result layer added on top of run_equating():
# conversions(), conversion_table(), equated(), and the summary()/plot() methods.

make_recipe <- function(seed = 42, design = "RG") {
  set.seed(seed)
  mk <- function(n, unique_items, common_items, p) {
    items <- c(unique_items, common_items)
    m <- matrix(stats::rbinom(n * length(items), 1, p), nrow = n,
                dimnames = list(NULL, items))
    cbind(data.frame(id = seq_len(n)), as.data.frame(m, check.names = FALSE))
  }
  common <- paste0("c", 1:5)
  fx <- mk(400, paste0("xA", 1:10), common, 0.55)
  fy <- mk(420, paste0("yB", 1:10), common, 0.50)

  init_equating() |>
    add_form(fx, name = "X", id_cols = "id", min_score = 0, max_score = 15, verbose = FALSE) |>
    add_form(fy, name = "Y", id_cols = "id", min_score = 0, max_score = 15, verbose = FALSE) |>
    add_plan(Y ~ X) |>
    add_design(design)
}

test_that("run_equating populates a tidy conversions table", {
  res <- make_recipe() |>
    add_method("identity") |>
    add_method("linear", boot_reps = 1) |>
    add_method("linear", mean_only = TRUE, boot_reps = 1) |>
    add_method("equipercentile") |>
    run_equating()

  conv <- conversions(res)
  expect_s3_class(conv, "data.frame")
  expect_true(all(c("from", "to", "method", "smooth", "x_score",
                    "equivalent_score", "se") %in% names(conv)))
  # Mean and (slope) Linear must both survive -- they used to collide on title.
  expect_true(all(c("Identity", "Linear", "Mean",
                    "Equipercentile (No Smoothing)") %in% conv$method))
  expect_setequal(unique(conv$from), "X")
  expect_setequal(unique(conv$to), "Y")
})

test_that("conversion_table() is wide with one column per method", {
  res <- make_recipe() |>
    add_method("identity") |>
    add_method("linear", boot_reps = 1) |>
    run_equating()

  ct <- conversion_table(res)
  expect_true("x_score" %in% names(ct))
  expect_true(all(c("Identity", "Linear") %in% names(ct)))
  expect_equal(ct$x_score, 0:15)
  # Identity maps a score to itself.
  expect_equal(ct$Identity, as.numeric(0:15))
  expect_identical(attr(ct, "from"), "X")
  expect_identical(attr(ct, "to"), "Y")
})

test_that("equated() applies a chosen conversion and validates the method", {
  res <- make_recipe() |>
    add_method("identity") |>
    add_method("linear", boot_reps = 1) |>
    run_equating()

  # At a tabulated score, equated() matches the conversion table exactly.
  ct <- conversion_table(res)
  expect_equal(equated(res, 10, method = "Linear"),
               ct$Linear[ct$x_score == 10])
  # Identity is the identity function.
  expect_equal(equated(res, c(3, 7, 12), method = "Identity"), c(3, 7, 12))
  # Interpolation between points stays within the neighbouring values.
  mid <- equated(res, 5.5, method = "Linear")
  expect_true(mid > ct$Linear[ct$x_score == 5] && mid < ct$Linear[ct$x_score == 6])
  # Unknown method is an error.
  expect_error(equated(res, 5, method = "Nope"), "not found")
})

test_that("summary() returns per-pairing tables and moments", {
  res <- make_recipe() |>
    add_method("linear", boot_reps = 1) |>
    add_method("equipercentile") |>
    run_equating()

  s <- summary(res)
  expect_type(s, "list")
  expect_true("X -> Y" %in% names(s))
  expect_s3_class(s[["X -> Y"]]$conversion_table, "data.frame")
  expect_true(is.matrix(s[["X -> Y"]]$moments))
})

test_that("a failing method does not abort the run", {
  # cubic-spline post-smoothing is unimplemented; it must fail in isolation.
  res <- make_recipe() |>
    add_method("equipercentile") |>
    add_method("equipercentile", smooth = "cubic_spline") |>
    run_equating()

  methods <- unique(conversions(res)$method)
  expect_true("Equipercentile (No Smoothing)" %in% methods)
})

test_that("common-item design yields Tucker/Levine/chained/FE conversions", {
  res <- make_recipe(design = "CG") |>
    add_method("linear", type = "all", boot_reps = 1) |>
    add_method("equipercentile", type = "E") |>
    run_equating()

  methods <- unique(conversions(res)$method)
  expect_true(all(c("Tucker", "Chained", "FrequencyEstimation") %in% methods))
  expect_length(equated(res, c(5, 10), method = "Tucker"), 2)
})

test_that("summary() includes a linear parameter comparison with SEs", {
  res <- make_recipe() |>
    add_method("linear", boot_reps = 50, boot_type = "perc") |>
    add_method("linear", mean_only = TRUE, boot_reps = 1) |>
    run_equating()

  s <- summary(res)
  params <- s[["X -> Y"]]$parameters
  expect_s3_class(params, "data.frame")
  expect_true(all(c("method", "slope", "se_slope", "intercept", "se_intercept") %in% names(params)))
  expect_true(all(c("Linear", "Mean") %in% params$method))
  # a bootstrapped linear method has a slope SE; mean equating fixes the slope to 1
  expect_false(is.na(params$se_slope[params$method == "Linear"]))
  expect_equal(params$slope[params$method == "Mean"], 1)
})

test_that("plot() renders default, difference, and SE-band views", {
  res <- make_recipe() |>
    add_method("linear", boot_reps = 50) |>
    add_method("identity") |>
    run_equating()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(res))
  expect_no_error(plot(res, difference = TRUE))
  expect_no_error(plot(res, difference = TRUE, se = TRUE))
})
