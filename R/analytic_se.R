# Analytic (delta-method) standard errors for equating.
# Port of the equipercentile case from Equating Recipes Analytic_SEs.c
# (Kolen & Brennan, Chapter 7 / Kolen & Brennan 2004, Section 7.3).

#' Analytic standard errors of equipercentile equating (random groups)
#'
#' Delta-method standard errors of the equipercentile Y-equivalents of the X
#' scores under the random-groups design. For score \eqn{x_i} with percentile
#' rank \eqn{P_i} (matched onto the Y scale), the standard error is
#' \deqn{\mathrm{se}[\hat e_Y(x_i)] = \frac{\sqrt{P_i(1-P_i)\,(1/N_X + 1/N_Y)}}{g_Y}}
#' where \eqn{g_Y} is the Y relative frequency (density) of the score interval
#' the rank maps into. A standard error of 0 is returned wherever the density is
#' 0 (e.g. zero-frequency scores), matching the original implementation.
#'
#' @param prdx Numeric vector of percentile ranks (0-100) for the X scores.
#' @param crfdy Cumulative relative frequency distribution of old form Y.
#' @param npx,npy Numbers of examinees for forms X and Y.
#' @return A numeric vector of estimated standard errors, one per X score.
#' @references Kolen, M. J., & Brennan, R. L. (2014). \emph{Test Equating,
#'   Scaling, and Linking} (3rd ed.), Section 7.3. Springer.
#' @export
se_ep_equate <- function(prdx, crfdy, npx, npy) {
  fdy <- diff(c(0, crfdy))                 # Y relative frequencies (density)
  P   <- prdx / 100
  yu  <- vapply(P, function(p) {
    if (!is.finite(p) || p <= 0) return(NA_integer_)
    w <- which(crfdy >= p - 1e-12)
    if (length(w)) w[1] else NA_integer_
  }, integer(1))
  gY <- ifelse(is.na(yu), 0, fdy[yu])
  ifelse(gY > 0 & P > 0 & P < 1, sqrt(P * (1 - P) * (1 / npx + 1 / npy)) / gY, 0)
}
