# Cubic smoothing spline postsmoothing (Kolen & Brennan "Equating Recipes",
# Chapter 11; Reinsch 1967). Seol is the author of the original C code.
#
# The smoothing spline g(x) minimizes the integrated squared second derivative
# subject to sum(((g(x_i) - y_i) / se_i)^2) <= S, fit over a restricted score
# range and linearly extended outside it. The scaled smoothing parameter `s`
# (Eq 11.4) maps to S = s * (number of nodes).

#' Fit a Reinsch smoothing cubic spline
#'
#' Solves for the smoothing cubic spline coefficients on `[x[1], x[N]]` using
#' the matrix formulation and Newton iteration for the Lagrange parameter from
#' Kolen & Brennan (Equating Recipes, Fig 11.1).
#'
#' @param x Strictly increasing node positions (length N).
#' @param y Values to smooth at the nodes.
#' @param se Standard errors (weights) at the nodes; larger se allows more
#'   deviation from `y`.
#' @param S Reinsch smoothing bound (`s * N` in the scaled parameterization).
#' @return A list of per-interval coefficients `a`, `b`, `c`, `d` (length N-1,
#'   with `a` the smoothed node values) plus the solved Lagrange parameter `p`.
#' @keywords internal
sspline <- function(x, y, se, S) {
  N <- length(x); n <- N - 1L; h <- diff(x); m <- N - 2L
  if (m < 1L) return(list(x = x, a = y, b = rep(0, max(n, 1)),
                          c = rep(0, max(n, 1)), d = rep(0, max(n, 1)), p = Inf))
  se[!is.finite(se) | se <= 0] <- min(se[is.finite(se) & se > 0], 1)
  D2 <- se^2

  Q  <- matrix(0, N, m)
  Tm <- matrix(0, m, m)
  for (j in seq_len(m)) {
    i <- j + 1L
    Q[i - 1L, j] <- 1 / h[i - 1L]
    Q[i,      j] <- -(1 / h[i - 1L] + 1 / h[i])
    Q[i + 1L, j] <- 1 / h[i]
    Tm[j, j] <- (2 / 3) * (h[i - 1L] + h[i])
    if (j < m) { Tm[j, j + 1L] <- h[i] / 3; Tm[j + 1L, j] <- h[i] / 3 }
  }
  QtD2Q <- crossprod(Q, D2 * Q)
  Qty   <- crossprod(Q, y)

  p <- 1; pnew <- 0; it <- 0L
  while (abs(pnew - p) > 1e-10 && it < 300L) {
    p <- pnew; it <- it + 1L
    R <- chol(QtD2Q + p * Tm)
    u <- backsolve(R, forwardsolve(t(R), Qty))
    e <- sum((se * (Q %*% u))^2)                 # F(p)^2 = ||D Q u||^2
    f <- sum(u * (Tm %*% u))                     # u' T u
    w <- forwardsolve(t(R), Tm %*% u)            # R' w = T u
    g <- sum(w^2)
    denom <- f - p * g
    if (!is.finite(denom) || abs(denom) < 1e-15) break
    pnew <- p + (e - sqrt(S * e)) / denom
    if (!is.finite(pnew) || pnew < 0) { pnew <- max(pnew, 0); break }
  }
  p <- max(pnew, 0)
  R <- chol(QtD2Q + p * Tm)
  u <- backsolve(R, forwardsolve(t(R), Qty))
  cc <- c(0, as.numeric(p * u), 0)               # c0 = cn = 0
  a  <- as.numeric(y - D2 * (Q %*% u))           # smoothed node values
  dd <- bb <- numeric(n)
  for (i in seq_len(n)) {
    dd[i] <- (cc[i + 1L] - cc[i]) / (3 * h[i])
    bb[i] <- (a[i + 1L] - a[i]) / h[i] - cc[i] * h[i] - dd[i] * h[i]^2
  }
  list(x = x, a = a, b = bb, c = cc[seq_len(n)], d = dd, p = p, a_all = a)
}

#' Cubic-spline postsmoothing of an equating function
#'
#' Fits a Reinsch smoothing spline (see [sspline()]) to the equipercentile
#' equivalents over the restricted range `xvalues[xlow:xhigh]` and evaluates it
#' at `vectX`, extending it with straight lines to the score-scale corners
#' `(-0.5, -0.5)` and `(KX + 0.5, KY + 0.5)` outside that range (Eq 11.5).
#'
#' @param xvalues Full vector of new-form (X) score points.
#' @param yvalues Unsmoothed equated (Y) equivalents at `xvalues`.
#' @param dyi Standard errors of the equated values at `xvalues`.
#' @param s Scaled smoothing parameter (Eq 11.4); `S = s * n_nodes`.
#' @param xlow,xhigh Integer index positions into `xvalues` bounding the fit
#'   (typically where the X percentile rank is within `[prlow, prhigh]`).
#' @param ky Maximum score on the Y scale (`KY`).
#' @param vectX Score points at which to return smoothed equivalents.
#' @return A list with `vectY` (smoothed equivalents at `vectX`) and the fitted
#'   spline `coefs`.
#' @keywords internal
post_smooth <- function(xvalues, yvalues, dyi, s, xlow, xhigh, ky, vectX) {
  idx <- xlow:xhigh
  xr <- xvalues[idx]; yr <- yvalues[idx]; sr <- dyi[idx]
  N  <- length(xr)
  if (N < 4L) {                                   # too few points to smooth
    return(list(vectY = stats::approx(xvalues, yvalues, xout = vectX, rule = 2)$y,
                coefs = NULL))
  }
  fit <- sspline(xr, yr, sr, s * N)

  KX <- max(xvalues); KY <- ky
  xlo <- xr[1]; xhi <- xr[N]; a1 <- fit$a[1]; aN <- fit$a[N]
  llow  <- (a1 + 0.5) / (xlo + 0.5)
  lhigh <- (aN - (KY + 0.5)) / (xhi - (KX + 0.5))

  vectY <- vapply(vectX, function(q) {
    if (q < xlo) return(llow * (q + 0.5) - 0.5)
    if (q > xhi) return(lhigh * (q - xhi) + aN)
    i <- max(which(xr <= q)); if (i >= N) i <- N - 1L
    dx <- q - xr[i]
    fit$a[i] + fit$b[i] * dx + fit$c[i] * dx^2 + fit$d[i] * dx^3
  }, numeric(1))

  list(vectY = vectY, coefs = fit)
}
