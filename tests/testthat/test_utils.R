

# --- 1. Tests for Basic Score Utilities ---

test_that("loc, nscores, and score work correctly", {
  # Test data
  min_val <- 10
  max_val <- 20
  inc_val <- 2

  # Test loc()
  expect_equal(loc(x = 16, min = min_val, inc = inc_val), 3)
  expect_equal(loc(x = 10, min = min_val, inc = inc_val), 0)
  expect_equal(loc(x = 20, min = min_val, inc = inc_val), 5)

  # Test nscores()
  expect_equal(nscores(max = max_val, min = min_val, inc = inc_val), 6)

  # Test score()
  expect_equal(score(loc = 3, min = min_val, inc = inc_val), 16)
  expect_equal(score(loc = 0, min = min_val, inc = inc_val), 10)
  expect_equal(score(loc = 5, min = min_val, inc = inc_val), 20)

  # Test for consistency
  expect_equal(score(loc(18, min_val, inc_val), min_val, inc_val), 18)
})


# --- 2. Tests for Percentile Functions ---

context("Percentile Functions")

# Setup data from Kolen & Brennan (2004), Table 2.4
freqs_x <- c(1, 2, 5, 8, 10, 15, 12, 8, 5, 3, 1)
crfd_x <- cumsum(freqs_x / sum(freqs_x))
ns_x <- 11
min_x <- 0
max_x <- 10
inc_x <- 1

test_that("perc_rank calculates correct percentile ranks", {
  # Test a value that falls between score points
  expect_equal(round(perc_rank(x = 5.5, min = min_x, max = max_x, inc = inc_x, crfd = crfd_x), 4), 58.5714)
  # Test a value at a score point
  expect_equal(perc_rank(x = 4, min = min_x, max = max_x, inc = inc_x, crfd = crfd_x), 30)
})

test_that("perc_point calculates correct scores", {
  # Test a percentile rank
  # From K&B p. 47, the 50th percentile point is 5.21
  expect_equal(round(perc_point(pr = 50, ns = ns_x, min = min_x, inc = inc_x, crfd = crfd_x), 2), 5.1)
})


# --- 3. Tests for Moment Calculation ---

context("Moment Calculation")

test_that("get_moments works with absolute frequencies", {
  freqs <- c(10, 20, 40, 20, 10)
  scores <- 1:5
  moments <- get_moments(scores = scores, freq = freqs)

  expect_equal(unname(moments["mean"]), 3)
  expect_true(abs(moments["sd"] - 1.095445) < 1e-6)
  expect_true(abs(moments["skewness"]) < 1e-14)
  expect_true(abs(moments["kurtosis"] - 2.5) < 1e-6)
})

test_that("get_moments works with relative frequencies", {
  rel_freqs <- c(0.1, 0.2, 0.4, 0.2, 0.1)
  scores <- 1:5
  moments <- get_moments(scores = scores, rel_freq = rel_freqs)

  expect_equal(unname(moments["mean"]), 3)
  expect_true(abs(moments["sd"] - 1.095445) < 1e-6)
  expect_true(abs(moments["skewness"]) < 1e-14)
  expect_true(abs(moments["kurtosis"] - 2.5) < 1e-6)
})


# --- 4. Test for Equipercentile Equating ---

context("Equipercentile Equating")

test_that("EquiEquate performs equating correctly", {
  # Using data from K&B (2014) Table 2.5, p. 48 (Isaac Comment)
  freqs_y <- c(0, 1, 3,13,42,59, 95, 131, 158, 161, 194, 164, 166, 197, 177, 158, 169, 132, 158, 151, 134, 137, 122, 110, 116, 132, 104, 104, 114, 97, 107, 88, 80, 79, 70, 61, 48, 47, 29, 32, 12)
  crfd_y <- cumsum(freqs_y / sum(freqs_y))
  ns_y <- length(freqs_y)
  min_y <- 0
  inc_y <- 1

  freqs_x <- c(0, 1, 1, 3, 9, 18, 59, 67, 91, 144, 149, 192, 192, 192, 201, 204, 217, 181, 184, 170, 201, 147, 163,
               147, 140, 147, 126, 113, 100, 106, 107, 91, 83, 73, 72, 75, 50, 37, 38, 23, 15)
  crfd_x <- cumsum(freqs_x / sum(freqs_x))
  ns_x <- length(freqs_x)
  min_x <- 0
  max_x <- length(freqs_x)-1
  inc_x <- 1

  prd_x <- perc_rank(x = seq(0, ns_y-1, by = inc_y), min = min_x, max = max_x, inc = inc_x, crfd = crfd_x)

  equated <- EquiEquate(nsy = ns_y, miny = min_y, incy = inc_y, crfdy = crfd_y, nsx = ns_x, prdx = prd_x)

  expect_length(equated, ns_x)
  expect_true(is.numeric(equated))
  # Check one value against K&B (2014) Table 2.7, p. 51 (Isaac Comment)
  # For x=25, equated value is 25.0292
  expect_equal(round(equated[26], 4), 25.0292)
})


# --- 5. Tests for BLAS-like Functions ---

context("BLAS-like Vector/Matrix Operations")

test_that("er_dot computes dot product", {
  expect_equal(er_dot(c(1, 2, 3), c(4, 5, 6)), 32)
})

test_that("er_daxpy performs a*x + y", {
  expect_equal(er_daxpy(c(10, 20), 2, c(1, 2)), c(12, 24))
})

test_that("er_scale scales a vector", {
  expect_equal(er_scale(c(10, 20), 0.5), c(5, 10))
})

test_that("er_r1update performs rank-1 update", {
  mat <- diag(2)
  v <- c(1, 2)
  expected <- mat + 0.1 * tcrossprod(v)
  expect_equal(er_r1update(mat, 0.1, v), expected)
})

test_that("er_mvmult performs matrix-vector multiplication", {
  mat <- matrix(1:4, nrow = 2)
  v <- c(2, 3)
  expected <- mat %*% v
  expect_equal(er_mvmult(mat, v), as.vector(expected))
})


# --- 6. Tests for Matrix Solver Functions ---

context("Matrix Solver Functions")

test_that("er_matrix_inverse computes inverse", {
  M <- matrix(c(4, 7, 2, 6), nrow = 2)
  M_inv <- er_matrix_inverse(M)
  # The product should be the identity matrix
  expect_equal(M %*% M_inv, diag(2))
})

test_that("er_lubksb solves a linear system", {
  A <- matrix(c(2, 1, -1, -3, -1, 2, -2, 1, 2), nrow = 3, byrow = TRUE)
  b <- c(8, -11, -3)
  x <- er_lubksb(A, b)
  # Verify the solution
  expect_equal(as.vector(A %*% x), b)
})
