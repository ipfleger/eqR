# Cubic-spline postsmoothing internals: analytic equipercentile SE and the
# Reinsch smoothing spline.

test_that("se_ep_equate reproduces K&B Table 7.2 (analytic equipercentile SE)", {
  Nx <- c(0,1,1,3,9,18,59,67,91,144,149,192,192,192,201,204,217,181,184,170,201,
          147,163,147,140,147,126,113,100,106,107,91,83,73,72,75,50,37,38,23,15)
  Ny <- c(0,1,3,13,42,59,95,131,158,161,194,164,166,197,177,158,169,132,158,151,
          134,137,122,110,116,132,104,104,114,97,107,88,80,79,70,61,48,47,29,32,12)
  NX <- sum(Nx); NY <- sum(Ny)
  Fx <- cumsum(Nx / NX); Gy <- cumsum(Ny / NY)
  prdx <- perc_rank(0:40, 0, 40, 1, Fx)
  se <- se_ep_equate(prdx, Gy, NX, NY)

  # K&B (2004) Table 7.2 analytic delta-method SEs.
  target <- c(0,0.83055,0.52100,0.82097,0.29502,0.14781,0.25411,0.15818,0.19691,
    0.17612,0.17312,0.19516,0.17995,0.23109,0.24312,0.21385,0.27635,0.26173,0.33835,
    0.28261,0.29473,0.32987,0.31827,0.38646,0.35546,0.30133,0.36831,0.35323,0.30691,
    0.34220,0.28963,0.32680,0.33093,0.30477,0.30798,0.30435,0.32400,0.27137,0.34301,
    0.20179,0.27872)

  expect_equal(se[1], 0)                                  # zero-frequency score
  # matches the published delta-method SEs across the interior score range
  mid <- 8:32
  expect_lt(max(abs(se[mid + 1] - target[mid + 1])), 0.02)
})

test_that("sspline preserves a straight line and smooths within range", {
  x <- 0:20
  # exactly linear data: a smoothing spline should reproduce it (zero roughness)
  y <- 2 + 0.5 * x
  out <- post_smooth(x, y, rep(1, length(x)), s = 0.2,
                     xlow = 1, xhigh = length(x), ky = 20, vectX = x)$vectY
  expect_equal(out, y, tolerance = 1e-3)

  # a spike gets pulled toward its neighbours (smoothed)
  y2 <- rep(5, length(x)); y2[11] <- 9
  out2 <- post_smooth(x, y2, rep(1, length(x)), s = 0.5,
                      xlow = 1, xhigh = length(x), ky = 20, vectX = x)$vectY
  expect_lt(out2[11], 9)
  expect_gt(out2[11], 5)
})
