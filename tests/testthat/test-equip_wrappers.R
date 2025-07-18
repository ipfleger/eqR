

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


context("Random Groups and Single Group Linear Equating")

test_that("linear_equate_rgsg calculates correctly for linear method", {
  mean_x <- 25
  sd_x <- 5
  mean_y <- 50
  sd_y <- 10

  result <- linear_equate_rgsg(
    mnx = mean_x, sdx = sd_x, mny = mean_y, sdy = sd_y,
    method = 'L',
    min_x = 10, max_x = 40, inc_x = 1
  )

  # Check slope and intercept
  expect_equal(result$a, 2.0) # 10 / 5
  expect_equal(result$b, 0.0) # 50 - 2 * 25

  # Check a few equated scores
  # For x = 10, y = 2*10 + 0 = 20
  # For x = 25, y = 2*25 + 0 = 50
  expect_equal(result$equated_scores[1], 20)
  expect_equal(result$equated_scores[16], 50)
})

test_that("linear_equate_rgsg calculates correctly for mean method", {
  mean_x <- 25
  sd_x <- 5 # Should be ignored
  mean_y <- 50
  sd_y <- 10 # Should be ignored

  result <- linear_equate_rgsg(
    mnx = mean_x, sdx = sd_x, mny = mean_y, sdy = sd_y,
    method = 'M',
    min_x = 10, max_x = 40, inc_x = 1
  )

  # Check slope and intercept
  expect_equal(result$a, 1.0) # Slope is fixed at 1 for mean equating
  expect_equal(result$b, 25.0) # 50 - 1 * 25

  # Check a few equated scores
  # For x = 10, y = 1*10 + 25 = 35
  # For x = 25, y = 1*25 + 25 = 50
  expect_equal(result$equated_scores[1], 35)
  expect_equal(result$equated_scores[16], 50)
})

test_that("linear_equate_rgsg handles errors", {
  # Test for invalid method
  expect_error(
    linear_equate_rgsg(25, 5, 50, 10, 'X', 10, 40, 1),
    "Invalid method specified"
  )

  # Test for zero standard deviation in linear equating
  expect_error(
    linear_equate_rgsg(25, 0, 50, 10, 'L', 10, 40, 1),
    "Standard deviation of X \\(sdx\\) cannot be zero"
  )
})

##### ------------


context("Random Groups and Single Group Linear Equating")

test_that("linear_equate_rgsg calculates correctly for linear method K&B 2.7", {

  moments_x <- get_moments(scores = 0:40, freq = c(0, 1, 1, 3, 9, 18, 59, 67, 91, 144, 149, 192, 192, 192, 201, 204, 217, 181, 184, 170, 201, 147, 163,
                                                   147, 140, 147, 126, 113, 100, 106, 107, 91, 83, 73, 72, 75, 50, 37, 38, 23, 15))

  moments_y <- get_moments(scores = 0:40, freq = c(0, 1, 3,13,42,59, 95, 131, 158, 161, 194, 164, 166, 197, 177, 158, 169, 132, 158, 151, 134, 137, 122, 110, 116, 132, 104, 104, 114, 97, 107, 88, 80, 79, 70, 61, 48, 47, 29, 32, 12))

  result <- linear_equate_rgsg(
    mnx = moments_x["mean"], sdx = moments_x["sd"], mny = moments_y["mean"], sdy = moments_y["sd"],
    method = 'L',
    min_x = 0, max_x = 40, inc_x = 1
  )

  m_result <- linear_equate_rgsg(
    mnx = moments_x["mean"], sdx = moments_x["sd"], mny = moments_y["mean"], sdy = moments_y["sd"],
    method = 'M',
    min_x = 0, max_x = 40, inc_x = 1
  )
  expect_equal(round(result$equated_scores, 4)[37], 36.5583)
  expect_equal(round(m_result$equated_scores, 4)[37], 35.1274)

})


test_that("linear_equate_rgsg calculates correctly for linear method", {
  mean_x <- 25
  sd_x <- 5
  mean_y <- 50
  sd_y <- 10

  result <- linear_equate_rgsg(
    mnx = mean_x, sdx = sd_x, mny = mean_y, sdy = sd_y,
    method = 'L',
    min_x = 10, max_x = 40, inc_x = 1
  )

  # Check slope and intercept
  expect_equal(result$a, 2.0) # 10 / 5
  expect_equal(result$b, 0.0) # 50 - 2 * 25

  # Check a few equated scores
  # For x = 10, y = 2*10 + 0 = 20
  # For x = 25, y = 2*25 + 0 = 50
  expect_equal(result$equated_scores[1], 20)
  expect_equal(result$equated_scores[16], 50)
})

test_that("linear_equate_rgsg calculates correctly for mean method", {
  mean_x <- 25
  sd_x <- 5 # Should be ignored
  mean_y <- 50
  sd_y <- 10 # Should be ignored

  result <- linear_equate_rgsg(
    mnx = mean_x, sdx = sd_x, mny = mean_y, sdy = sd_y,
    method = 'M',
    min_x = 10, max_x = 40, inc_x = 1
  )

  # Check slope and intercept
  expect_equal(result$a, 1.0) # Slope is fixed at 1 for mean equating
  expect_equal(result$b, 25.0) # 50 - 1 * 25

  # Check a few equated scores
  # For x = 10, y = 1*10 + 25 = 35
  # For x = 25, y = 1*25 + 25 = 50
  expect_equal(result$equated_scores[1], 35)
  expect_equal(result$equated_scores[16], 50)
})

test_that("linear_equate_rgsg handles errors", {
  # Test for invalid method
  expect_error(
    linear_equate_rgsg(25, 5, 50, 10, 'X', 10, 40, 1),
    "Invalid method specified"
  )

  # Test for zero standard deviation in linear equating
  expect_error(
    linear_equate_rgsg(25, 0, 50, 10, 'L', 10, 40, 1),
    "Standard deviation of X \\(sdx\\) cannot be zero"
  )
})

# _------------- CINEG Linear equating


context("CINEG Linear Equating")

# Mock data based on Kolen & Brennan (2004), Table 4.1, p. 100
# Population 1 (Form X)
mnx1 <- 15.8205; sdx1 <- 6.5278; mnv1 <- 5.1063; sdv1 <- 2.3760
# Population 2 (Form Y)
mny2 <- 18.6728; sdy2 <- 6.8784; mnv2 <- 5.8626; sdv2 <- 2.4515
# Covariances (calculated from correlation r=0.80)
covxv1 <- 13.4088 # 0.8645 * sdx1 * sdv1
covyv2 <- 14.7603 # 0.8753 * sdy2 * sdv2
# Other parameters
w1 <- 1; anchor <- 1 # External anchor
min_x <- 0; max_x <- 36; inc_x <- 1

test_that("linear_equate_ci calculates all methods correctly", {

  results <- linear_equate_ci(
    mnx1, sdx1, mnv1, sdv1, covxv1,
    mny2, sdy2, mnv2, sdv2, covyv2,
    w1, anchor, 'L', min_x, max_x, inc_x
  )

  # --- Check against K&B (2014) Eq. 4.82, p. 125 ---

  # Tucker
  expect_equal(round(results$summary$a[1], 4), 1.0289)
  expect_equal(round(results$summary$b[1], 4), 0.5370)

  # --- Check against K&B (2014) Eq. 4.83, p. 126 ---

  # Levine Observed
  expect_equal(round(results$summary$a[2], 3), 1.011) # I think there is some rounding error here as opposed to calculation error
  expect_equal(round(results$summary$b[2], 4), 0.2517)

  # Levine True
  expect_equal(round(results$summary$a[3], 3), 1.009)
  expect_equal(round(results$summary$b[3], 3), 0.291)

  # Chained
  expect_equal(round(results$summary$a[4], 4), 1.0213)
  expect_equal(round(results$summary$b[4], 4), 0.3940)

  # Check an equated score
  # For Tucker, equated score for x=20 is 1.0289*20 + 0.5370 = 21.115
  expect_equal(round(results$equated_scores[1, "x_20"], 2), 21.12)
  expect_equal(round(results$equated_scores[4, "x_20"], 2), 20.82)
  expect_equal(round(results$equated_scores[2, "x_20"], 2), 20.47)
  expect_equal(round(results$equated_scores[3, "x_20"], 2), 20.46)
})

#---------------- Kernel Equating


context("Kernel Equating Functions")

# --- Mock Data ---
# A simple, symmetric distribution for testing
scores <- 0:10
freqs <- c(1, 2, 5, 8, 10, 15, 10, 8, 5, 2, 1)
rel_freqs <- freqs / sum(freqs)
hx <- 1.5 # A plausible bandwidth

test_that("std_normal_pdf and std_normal_cdf match base R", {
  expect_equal(std_normal_pdf(0), dnorm(0))
  expect_equal(std_normal_pdf(1.96), dnorm(1.96))

  expect_equal(std_normal_cdf(0), pnorm(0))
  expect_equal(std_normal_cdf(1.96), pnorm(1.96))
})

test_that("kernel_continu_pdf runs and returns a single numeric value", {
  pdf_val <- kernel_continu_pdf(x = 5, scores = scores, rel_freq = rel_freqs, hx = hx)
  expect_true(is.numeric(pdf_val))
  expect_length(pdf_val, 1)
  # The peak of the PDF should be near the mean
  mean_val <- sum(scores * rel_freqs)
  expect_gt(
    kernel_continu_pdf(mean_val, scores, rel_freqs, hx),
    kernel_continu_pdf(mean_val + 3, scores, rel_freqs, hx)
  )
})

test_that("kernel_continu_cdf runs and returns a valid probability", {
  cdf_val <- kernel_continu_cdf(x = 5, scores = scores, rel_freq = rel_freqs, hx = hx)
  expect_true(is.numeric(cdf_val))
  expect_length(cdf_val, 1)
  expect_gte(cdf_val, 0)
  expect_lte(cdf_val, 1)

  # CDF should be monotonic
  expect_lt(
    kernel_continu_cdf(4, scores, rel_freqs, hx),
    kernel_continu_cdf(6, scores, rel_freqs, hx)
  )
})

test_that("kernel_inverse_cdf finds the correct score for a given probability", {
  # The CDF at the mean should be around 0.5 for a symmetric distribution
  mean_val <- sum(scores * rel_freqs)
  cdf_at_mean <- kernel_continu_cdf(mean_val, scores, rel_freqs, hx)

  # Find the score that corresponds to the CDF value at the mean
  inverted_score <- kernel_inverse_cdf(p = cdf_at_mean, scores = scores, rel_freq = rel_freqs, hx = hx)

  # It should be very close to the actual mean
  expect_equal(inverted_score, mean_val, tolerance = 1e-6)

  # Test the median (50th percentile)
  median_score <- kernel_inverse_cdf(p = 0.5, scores = scores, rel_freq = rel_freqs, hx = hx)
  expect_equal(median_score, mean_val, tolerance = 1e-6)
})


#-----------LogLinear--------------

context("Log-Linear Smoothing Functions")

# --- Mock Data ---
# A simple univariate distribution
freqs <- c(10, 20, 40, 20, 10)
scores <- 1:5
n_total <- sum(freqs)
n_scores <- length(scores)

# --- Tests for Helper Functions ---

test_that("design_matrix creates a matrix of correct dimensions", {
  design <- design_matrix(nsu = n_scores, minu = 1, incu = 1, cu = 3)
  expect_equal(dim(design$B), c(n_scores, 3))
  expect_equal(design$B_raw[, 1], scores)
  expect_equal(design$B_raw[, 2], scores^2)
})

test_that("get_LLmoments calculates moments correctly", {
  design <- design_matrix(nsu = n_scores, minu = 1, incu = 1, cu = 2)
  moments <- get_LLmoments(design$B, design$B_raw, freqs, n_total)

  # The first moment should be the mean
  expect_equal(moments$mts_raw[1], 3)
})

# --- Test for the Main Iteration Function ---

test_that("iteration function runs and converges", {
  design <- design_matrix(nsu = n_scores, minu = 1, incu = 1, cu = 2, scale = TRUE)

  # Expect the iteration to run without error
  results <- expect_no_error(
    iteration(B = design$B, B_raw = design$B_raw, nct = freqs, N = n_total)
  )

  # Check that the results have the expected structure
  expect_true(is.list(results))
  expect_named(results, c("Beta", "mct", "nit", "lrchisq", "n_mts", "m_mts", "ap", "converged"))

  # The sum of fitted frequencies should equal the total N
  expect_equal(sum(results$mct), n_total, tolerance = 1e-6)

  # The fitted moments should be very close to the observed moments
  expect_equal(results$n_mts, results$m_mts, tolerance = 1e-4)
})

test_that("smooth_ull wrapper runs correctly", {

  # Expect the full smoothing wrapper to run without error
  results <- expect_no_error(
    smooth_ull(n = n_total, ns = n_scores, min = 1, inc = 1, fd = freqs, c = 2)
  )

  # Check for the final outputs
  expect_true(is.list(results))
  expect_true("density" %in% names(results))
  expect_true("crfd" %in% names(results))
  expect_true("prd" %in% names(results))
  expect_equal(sum(results$density), 1.0, tolerance = 1e-6)
})

# ------------------ CG_EquiEquate ---------------

context("CINEG Equating with Textbook Data")

# --- Data from Kolen & Brennan (2014), Chapter 5 ---

# Group A (took Form A and Anchor V)
bxvin <- rbind(
  c(.04, .04, .02, .00),
  c(.04, .08, .02, .01),
  c(.06, .12, .05, .02),
  c(.03, .12, .05, .05),
  c(.02, .03, .04, .06),
  c(.01, .01, .02, .06)
)

# Group B (took Form B and Anchor V)
byvin <- rbind(
  c(.04, .03, .01, .00),
  c(.07, .05, .07, .01),
  c(.03, .05, .12, .02),
  c(.03, .04, .13, .05),
  c(.02, .02, .05, .06),
  c(.01, .01, .02, .06)
)

# Other parameters from the example
w1 <- 1
internal <- FALSE
rv1 <- .8645; rv2 <- .8753

# Score scales
minx <- 0; maxx <- 5; incx <- 1; nsx <- 6
miny <- 0; maxy <- 5; incy <- 1; nsy <- 6
minv <- 0; maxv <- 3; incv <- 1; nsv <- 4


test_that("fe_mfe_equate produces correct FE and BH results", {

  # This test uses rv1=0, rv2=0 to specify the Frequency Estimation (FE) method.
  freq_results <- fe_mfe_equate(
    w1 = w1, internal = internal, bxvin = bxvin, byvin = byvin,
    rv1 = 0, rv2 = 0,
    minx = minx, maxx = maxx, incx = incx,
    miny = miny, maxy = maxy, incy = incy
  )

  # We need to specify the correlations to get the Braun Holland Method
  results <- fe_mfe_equate(
    w1 = w1, internal = internal, bxvin = bxvin, byvin = byvin,
    rv1 = rv1, rv2 = rv2,
    minx = minx, maxx = maxx, incx = incx,
    miny = miny, maxy = maxy, incy = incy
  )

  # Expected results from K&B (2014), p. 149
  expected_eraw <- c(-0.02, 0.83, 1.76, 2.92, 3.98, 5.0)
  expected_a_bh <- 1.0331
  expected_b_bh <- -0.1927

  # Check equipercentile results
  expect_equal(freq_results$eraw, expected_eraw, tolerance = 1e-2)

  # Check Braun-Holland parameters
  expect_equal(results$eraw_bh, c(-.1927, .8404, 1.8735, 2.9066, 3.9397, 4.9728), tolerance = 1e-4)
  expect_equal(results$a_bh, expected_a_bh, tolerance = 1e-4)
  expect_equal(results$b_bh, expected_b_bh, tolerance = 1e-4)
})


test_that("chained_equate with smoothing produces correct intermediate and final results", {

  # Calculate necessary inputs from the provided distributions
  marginal_a_x <- rowSums(bxvin)
  marginal_a_v <- colSums(bxvin)
  marginal_b_y <- rowSums(byvin)
  marginal_b_v <- colSums(byvin)

  crfd_a_x <- cumsum(marginal_a_x)
  crfd_a_v <- cumsum(marginal_a_v)
  crfd_b_y <- cumsum(marginal_b_y)
  crfd_b_v <- cumsum(marginal_b_v)

  # Get percentile ranks of Form A scores in Group A
  prdx1 <- perc_rank(x = 0:5, min = minx, max = maxx, inc = incx, crfd = crfd_a_x)

  # --- Test the intermediate step: Equating X to V ---
  extov <- EquiEquate(nsy = nsv, miny = minv, incy = incv, crfdy = crfd_a_v, nsx = nsx, prdx = prdx1)
  expected_extov <- c(-0.2500, 0.3750, 0.9375, 1.6250, 2.6250, 3.2500)
  expect_equal(extov, expected_extov, tolerance = 1e-4)

  # --- Test the final chained equating result ---
  final_equated <- chained_equate(
    nsx = nsx, prx1 = prdx1,
    minv = minv, maxv = maxv, incv = incv, nsv = nsv, Gv1 = crfd_a_v,
    miny = miny, incy = incy, nsy = nsy, Gy2 = crfd_b_y, Gv2 = crfd_b_v
  )

  expect_length(final_equated, nsx)
  expect_true(is.numeric(final_equated))

  # Based on manual calculation, the first value should be 0.125
  expect_equal(final_equated[1], 0.125)
})

# ---------- Bootstrap --------------------------------


context("Bootstrap Functions")
# I didn't check these for accuracy, I just wanted these to run correctly
# --- Mock Data ---
set.seed(456)
n_total <- 1000
original_freqs <- c(50, 150, 300, 300, 150, 50)
original_bfd <- matrix(c(10,20,30, 40,100,150, 60,150,100, 40,100,150, 30,20,10), nrow=5)
n_bfd <- sum(original_bfd)

test_that("bootstrap_univariate returns a valid frequency distribution", {
  boot_freq <- bootstrap_univariate(n = n_total, freq = original_freqs)

  expect_length(boot_freq, length(original_freqs))
  expect_equal(sum(boot_freq), n_total)
  expect_true(all(boot_freq >= 0))
})

test_that("bootstrap_bivariate returns a valid bivariate distribution", {
  boot_bfd <- bootstrap_bivariate(n = n_bfd, bfd = original_bfd)

  expect_equal(dim(boot_bfd), dim(original_bfd))
  expect_equal(sum(boot_bfd), n_bfd)
  expect_true(all(boot_bfd >= 0))
})

test_that("bootstrap_parametric_univ returns a valid frequency distribution", {
  smoothed_probs <- original_freqs / n_total
  boot_freq <- bootstrap_parametric_univ(n = n_total, smoothed_dist = smoothed_probs)

  expect_length(boot_freq, length(smoothed_probs))
  expect_equal(sum(boot_freq), n_total)
})

test_that("calculate_bootstrap_se computes summary statistics correctly", {
  # Create a matrix of mock bootstrap results
  # 100 replications for a test with 5 score points
  results <- matrix(rnorm(5 * 100, mean = 10, sd = 2), nrow = 5, ncol = 100)

  summary_stats <- calculate_bootstrap_se(results, original_freq = c(10,20,40,20,10))

  expect_length(summary_stats$mean_scores, 5)
  expect_length(summary_stats$se_scores, 5)
  expect_true(is.numeric(summary_stats$overall_se))

  # The mean of the bootstrap means should be close to the true mean
  expect_equal(mean(summary_stats$mean_scores), 10, tolerance = 0.5)
  # The mean of the bootstrap SEs should be close to the true SD
  expect_equal(mean(summary_stats$se_scores), 2, tolerance = 0.5)
})


# ---------- Advanved Utility Functions ---------------


context("Advanced Utility Functions")

test_that("cum_rel_freqs works correctly", {
  rel_freqs <- c(0.1, 0.2, 0.3, 0.4)
  expected <- c(0.1, 0.3, 0.6, 1.0)
  expect_equal(cum_rel_freqs(rel_freqs), expected)
})

test_that("equated_ss performs interpolation and rounding", {
  # Raw-to-scale conversion table for Y
  y_conversion_table <- data.frame(
    raw = c(0, 10, 20, 30),
    scale = c(100, 200, 300, 400)
  )

  # Equated raw scores on Y's scale
  equated_raw_scores <- c(5, 15, 25)

  results <- equated_ss(
    eraw = equated_raw_scores,
    yct = y_conversion_table,
    lprss = 100,
    hprss = 400,
    round_digits = 0
  )

  # Check unrounded scores (should be 150, 250, 350)
  expect_equal(results$essu, c(150, 250, 350))
  # Check rounded scores
  expect_equal(results$essr, c(150, 250, 350))
})

test_that("er_qrdcmp performs QR decomposition", {
  mat <- matrix(c(1, -1, 4, 1, 4, -2, 1, 4, 2, 1, -1, 0), nrow = 4)
  qr_result <- er_qrdcmp(mat)

  expect_s3_class(qr_result, "qr")
  # Reconstruct the original matrix to verify
  Q <- qr.Q(qr_result)
  R <- qr.R(qr_result)
  expect_equal(Q %*% R, mat, tolerance = 1e-9)
})

test_that("er_rtsafe finds the root of a function", {
  # Find the root of cos(x) in the interval [0, 3]
  root <- er_rtsafe(cos, interval = c(0, 4))
  expect_equal(root, pi / 2, tolerance = 1e-7)
})

test_that("er_dfpmin minimizes a function", {
  # Rosenbrock banana function
  rosenbrock <- function(x) {
    (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  }

  # Start optimization from a point
  result <- er_dfpmin(par = c(-1.2, 1), fn = rosenbrock)

  # The minimum should be at (1, 1) with a value of 0
  expect_equal(result$value, 0, tolerance = 1e-6)
  expect_equal(result$par, c(1, 1), tolerance = 1e-4)
})


# You would need to install the VGAM package first:
# install.packages("VGAM")

# Assume our translated obs_density is in the environment

test_that("obs_density matches a known beta-binomial implementation", {

  # Define parameters
  n_items <- 20
  alpha <- 5
  beta_p <- 3

  # Our function's output (as probabilities)
  beta_params <- c(alpha = alpha, beta = beta_p, lower = 0, upper = 1)
  our_counts <- obs_density(n_items, n_persons = 1, beta_params) # n_persons=1 for probabilities
  our_probs <- our_counts / sum(our_counts)
  plot(our_probs)
  # The trusted package's output
  # Note: VGAM's function takes size (n_items), shape1 (alpha), and shape2 (beta)
  vgam_probs <- VGAM::dbetabinom.ab(x = 0:n_items, size = n_items, shape1 = alpha, shape2 = beta_p)

  # Compare the two distributions
  expect_equal(our_probs, vgam_probs, tolerance = 1e-6)
})



# ------------- CLL -----------------
# Title: Bivariate Frequency Matrix from Equating Recipes Manual
#
# Description:
# This script creates an R matrix containing the bivariate frequency
# distribution of scores on the non-common items (X') versus scores on
# the common items (V), as shown in Table 10.14 of the Equating
# Recipes manual.

# Create the matrix by entering the data row by row
bivariate_freq_table <- matrix(
  c(
    # V=0 to V=12
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # X'=0
    0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # X'=1
    1, 2, 3, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0,  # X'=2
    2, 4, 5, 15, 6, 3, 2, 0, 0, 0, 0, 0, 0, # X'=3
    0, 4, 13, 11, 16, 2, 2, 2, 1, 0, 0, 0, 0, # X'=4
    2, 7, 19, 30, 18, 9, 9, 1, 0, 1, 0, 0, 0, # X'=5
    0, 9, 19, 29, 30, 20, 16, 0, 0, 0, 0, 0, 0, # X'=6
    6, 7, 20, 30, 33, 29, 11, 5, 2, 0, 1, 0, 0, # X'=7
    0, 12, 16, 29, 36, 24, 19, 8, 1, 0, 0, 0, 0, # X'=8
    3, 5, 18, 29, 28, 27, 15, 10, 3, 0, 0, 0, 0, # X'=9
    0, 2, 16, 21, 31, 23, 30, 12, 6, 1, 0, 0, 0, # X'=10
    0, 0, 4, 26, 20, 25, 24, 11, 6, 2, 1, 0, 0, # X'=11
    0, 2, 4, 11, 16, 27, 26, 9, 10, 2, 1, 0, 0, # X'=12
    0, 0, 1, 7, 16, 16, 17, 24, 8, 3, 2, 0, 0, # X'=13
    0, 0, 2, 3, 10, 12, 16, 20, 11, 6, 1, 0, 0, # X'=14
    0, 0, 1, 4, 5, 16, 13, 20, 7, 7, 2, 0, 0, # X'=15
    0, 0, 0, 0, 5, 7, 16, 17, 13, 9, 1, 3, 0, # X'=16
    0, 0, 0, 0, 3, 3, 7, 12, 16, 9, 3, 0, 0, # X'=17
    0, 0, 0, 0, 0, 2, 2, 11, 16, 14, 3, 4, 1, # X'=18
    0, 0, 0, 0, 0, 0, 2, 5, 7, 12, 4, 2, 0, # X'=19
    0, 0, 0, 0, 0, 2, 0, 1, 6, 7, 12, 3, 1, # X'=20
    0, 0, 0, 0, 0, 0, 3, 2, 4, 1, 4, 8, 2, # X'=21
    0, 0, 0, 0, 0, 0, 0, 3, 1, 1, 5, 3, 3, # X'=22
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 0, # X'=23
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1  # X'=24
  ),
  nrow = 25,
  ncol = 13,
  byrow = TRUE
)

# Add row and column names for clarity
rownames(bivariate_freq_table) <- paste0("X'=", 0:24)
colnames(bivariate_freq_table) <- paste0("V=", 0:12)

# Print the matrix to verify
print(bivariate_freq_table)

#         V_Score 0 V_Score 1 V_Score 2 V_Score 3
# A_Score 0        40        40        20         0
# A_Score 1        40        80        20        10
# A_Score 2        60       120        50        20
# A_Score 3        30       120        50        50
# A_Score 4        20        30        40        60
# A_Score 5        10        10        20        60

# The sum of all cells is the total number of people
print(sum(bivariate_freq_table))
# [1] 1000


n <- sum(bivariate_freq_table)
nsu <- nrow(bivariate_freq_table)
min_u <- 0
incu <- 1
nsv <- ncol(bivariate_freq_table)
minv <- 0
incv <- 1
nct <- as.vector((bivariate_freq_table))
scale = TRUE# Set to FALSE for continuized_log_linear, TRUE for all others


test <- smooth_bll(n = n, nsu = nsu, minu = min_u, incu = incu, nsv = nsv,
                   minv = minv, incv = incv, nct = nct, cu = 6, cv = 6,
                   cuv = 1, cpm = matrix(c(1,1), nrow = 1), scale = scale, crit = 1e-05)


test_output <- matrix(test$mct,
                      nrow = nrow(bivariate_freq_table),
                      ncol = ncol(bivariate_freq_table),
                      byrow = FALSE) |>
  round(3)

# Title: Fitted Bivariate Frequency Matrix from Equating Recipes Manual
#
# Description:
# This script creates an R matrix containing the log-linear smoothed
# (fitted) bivariate frequency distribution of scores on the non-common
# items (X') versus scores on the common items (V), as shown in
# Table 10.15 of the Equating Recipes manual.

# Create the matrix by entering the data row by row
fitted_freq_table <- matrix(
  c(
    # V=0 to V=12
    0.041, 0.092, 0.106, 0.070, 0.030, 0.009, 0.002, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.243, 0.618, 0.803, 0.603, 0.291, 0.097, 0.024, 0.004, 0.001, 0.000, 0.000, 0.000, 0.000,
    0.779, 2.236, 3.285, 2.791, 1.521, 0.576, 0.160, 0.033, 0.005, 0.001, 0.000, 0.000, 0.000,
    1.559, 5.056, 8.397, 8.067, 4.971, 2.129, 0.668, 0.158, 0.028, 0.004, 0.000, 0.000, 0.000,
    2.180, 7.993, 15.010, 16.302, 11.358, 5.499, 1.951, 0.521, 0.106, 0.016, 0.002, 0.000, 0.000,
    2.322, 9.627, 20.438, 25.097, 19.770, 10.820, 4.340, 1.312, 0.300, 0.052, 0.007, 0.001, 0.000,
    2.009, 9.416, 22.603, 31.379, 27.947, 17.294, 7.843, 2.680, 0.694, 0.135, 0.020, 0.002, 0.000,
    1.479, 7.840, 21.278, 33.397, 33.629, 23.527, 12.063, 4.660, 1.365, 0.300, 0.049, 0.006, 0.001,
    0.959, 5.748, 17.635, 31.295, 35.628, 28.181, 16.336, 7.135, 2.362, 0.588, 0.109, 0.015, 0.002,
    0.561, 3.800, 13.181, 26.445, 34.038, 30.440, 19.950, 9.851, 3.687, 1.038, 0.217, 0.034, 0.004,
    0.301, 2.304, 9.035, 20.495, 29.825, 30.155, 22.345, 12.475, 5.279, 1.679, 0.397, 0.069, 0.009,
    0.150, 1.297, 5.749, 14.745, 24.259, 27.731, 23.232, 14.664, 7.016, 2.524, 0.674, 0.133, 0.020,
    0.070, 0.683, 3.426, 9.935, 18.481, 23.885, 22.623, 16.145, 8.734, 3.551, 1.073, 0.240, 0.041,
    0.031, 0.340, 1.925, 6.312, 13.275, 19.397, 20.772, 16.760, 10.250, 4.713, 1.609, 0.407, 0.078,
    0.013, 0.160, 1.025, 3.800, 9.036, 14.928, 18.074, 16.487, 11.400, 5.926, 2.288, 0.655, 0.142,
    0.005, 0.072, 0.519, 2.175, 5.847, 10.922, 14.950, 15.419, 12.054, 7.084, 3.092, 1.000, 0.246,
    0.002, 0.031, 0.250, 1.184, 3.600, 7.601, 11.764, 13.718, 12.125, 8.056, 3.975, 1.454, 0.404,
    0.001, 0.012, 0.114, 0.612, 2.101, 5.017, 8.779, 11.573, 11.565, 8.688, 4.847, 2.004, 0.630,
    0.000, 0.005, 0.049, 0.297, 1.153, 3.113, 6.158, 9.178, 10.369, 8.807, 5.555, 2.597, 0.923,
    0.000, 0.002, 0.019, 0.133, 0.585, 1.785, 3.992, 6.726, 8.592, 8.251, 5.884, 3.110, 1.249,
    0.000, 0.001, 0.007, 0.054, 0.266, 0.919, 2.325, 4.430, 6.397, 6.945, 5.600, 3.346, 1.520,
    0.000, 0.000, 0.002, 0.019, 0.105, 0.408, 1.166, 2.511, 4.100, 5.033, 4.588, 3.100, 1.591,
    0.000, 0.000, 0.001, 0.005, 0.033, 0.146, 0.473, 1.152, 2.127, 2.951, 3.042, 2.323, 1.349,
    0.000, 0.000, 0.000, 0.001, 0.008, 0.039, 0.143, 0.393, 0.819, 1.285, 1.498, 1.294, 0.849,
    0.000, 0.000, 0.000, 0.000, 0.001, 0.007, 0.028, 0.089, 0.209, 0.371, 0.489, 0.477, 0.354
  ),
  nrow = 25,
  ncol = 13,
  byrow = TRUE
)

# Add row and column names for clarity
rownames(fitted_freq_table) <- paste0("X'=", 0:24)
colnames(fitted_freq_table) <- paste0("V=", 0:12)

# Print the matrix to verify
print(fitted_freq_table)

expect_true(all(abs(test_output - fitted_freq_table) < 3e-3))



# 1. Start with your raw data: a bivariate frequency table
# (e.g., scores on Form X vs. scores on Form Y from the same group)
bivariate_counts <- matrix(bivariate_freq_table)
total_n <- sum(bivariate_counts)

# 2. Call smooth_bll() to perform the log-linear smoothing.
#    This fits the model and returns a complex list object.
bivar_list <- smooth_bll(
  n = total_n,
  nsu = nrow(bivariate_counts), # Number of X scores
  nsv = ncol(bivariate_counts), # Number of Y scores
  nct = as.vector(t(bivariate_counts)), # The flattened frequency table
  cu = 3
)

# 3. Pass the resulting 'bivar_list' object to equate_cll().
#    The equate_cll function uses the parameters inside bivar_list
#    to perform the final equating calculations.
equated_scores <- equate_cll(design = "SG", bivar = bivar_list)

# The 'equated_scores' vector now holds your final result.
print(equated_scores)
#-----------------------------

# --- 1. SETUP ---
# Load necessary libraries and source the translated R functions.
# Make sure these files are in your working directory.
library(MASS)

# --- 2. LOAD DATA ---
# This is the ACT Math data from "actmathfreq.dat", used in the manual's example.
# We'll create the data frames here for a reproducible example.
act_data <- data.frame(
  score = 0:40,
  freq_x = c(0, 1, 1, 3, 9, 18, 59, 67, 91, 144, 149, 192, 192, 192, 201, 204, 217, 181, 184, 170, 201, 147, 163, 147, 140, 147, 126, 113, 100, 106, 107, 91, 83, 73, 72, 75, 50, 37, 38, 23, 15),
  freq_y = c(0, 1, 3, 13, 42, 59, 95, 131, 158, 161, 194, 164, 166, 197, 177, 158, 169, 132, 158, 151, 134, 137, 122, 110, 116, 132, 104, 104, 114, 97, 107, 88, 80, 79, 70, 61, 48, 47, 29, 32, 12)
)

n_x <- sum(act_data$freq_x)
n_y <- sum(act_data$freq_y)

# --- 3. LOG-LINEAR PRE-SMOOTHING ---
# As per the manual, we first smooth each distribution with a 6th-degree polynomial.
cat("--- Performing Log-Linear Smoothing ---\n")
smooth_x <- smooth_ull(n = n_x, ns = 41, min = 0, inc = 1, fd = act_data$freq_x, c = 6, scale = TRUE)
smooth_y <- smooth_ull(n = n_y, ns = 41, min = 0, inc = 1, fd = act_data$freq_y, c = 6, scale = TRUE)

# The kernel method uses the smoothed relative frequencies (density)
rel_freq_x_smoothed <- smooth_x$density
rel_freq_y_smoothed <- smooth_y$density

# --- 4. FIND OPTIMAL BANDWIDTH (h) ---
cat("--- Finding Optimal Bandwidth for Form X ---\n")
hx <- find_optimal_h(scores = act_data$score, rel_freq = rel_freq_x_smoothed)
cat("Optimal h for X:", hx, "\n")

cat("--- Finding Optimal Bandwidth for Form Y ---\n")
hy <- find_optimal_h(scores = act_data$score, rel_freq = rel_freq_y_smoothed)
cat("Optimal h for Y:", hy, "\n")

# --- 5. PERFORM KERNEL EQUATING ---
cat("--- Performing Kernel Equating ---\n")
equated_scores <- kernel_equate(
  scores_x = act_data$score,
  rel_freq_x = rel_freq_x_smoothed,
  hx = hx,
  scores_y = act_data$score,
  rel_freq_y = rel_freq_y_smoothed,
  hy = hy
)

# --- 6. VALIDATE RESULTS ---
results_df <- data.frame(
  x_score = act_data$score,
  equated_y_score = equated_scores
)

P <- act_data$freq_x
Q <- act_data$freq_y
xscores <- act_data$score

# 2. Calculate Probabilities and Sample Sizes
# The sample size for each group (N for test X, M for test Y)
N <- sum(P)
M <- sum(Q)

act_equating_smooth <- kequate::kequate(design = "EG",
                                        x = xscores,
                                        y = xscores,
                                        r = rel_freq_x_smoothed,
                                        s = rel_freq_y_smoothed,
                                        N = N,
                                        M = M,
                                        smoothed = FALSE)

# Ok, so we are
all(act_equating_smooth@equating$eqYx - equated_scores < 1e-4)

cat("\n--- Equating Results ---\n")
print(head(results_df))
cat("...\n")
print(tail(results_df))

cat("\nVALIDATION: Compare the 'equated_y_score' column with the results in Table 12.2 (page 181) of ER.pdf.\n")
# For example, the value for x_score = 20 should be close to 19.24488

# --- 1. Define a function to calculate raw moments ---
# This function calculates the first 'k' raw moments of a distribution
get_raw_moments <- function(scores, rel_freq, k = 6) {
  sapply(1:k, function(power) {
    sum(rel_freq * (scores^power))
  })
}

# --- 2. Get moments from the ORIGINAL data ---
original_rel_freq_x <- act_data$freq_x / sum(act_data$freq_x)
original_moments <- get_raw_moments(act_data$score, original_rel_freq_x)

# --- 3. Get moments from our SMOOTHED data ---
# (Using the smooth_x object from our previous test script)
smoothed_moments <- get_raw_moments(act_data$score, smooth_x$density)

# --- 4. Compare the moments ---
moment_comparison <- data.frame(
  Moment = 1:6,
  Original_Data = original_moments,
  Smoothed_Data = smoothed_moments,
  Difference = original_moments - smoothed_moments
)

cat("--- Comparison of the First 6 Raw Moments ---\n")
print(moment_comparison, digits = 10)
