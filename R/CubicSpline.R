#' Perform Cholesky Decomposition
#'
#' @description
#' This function computes the Cholesky decomposition of a symmetric,
#' positive-definite matrix. It serves as an R equivalent to the `dpdch`
#' function from the "Equating Recipes" C library.
#'
#' @details
#' The function uses R's built-in `chol()` function, which is highly optimized
#' for this purpose. The original C function `dpdch` overwrites the input
#' matrix, whereas this R function returns a new matrix containing the upper
#' triangular factor of the decomposition.
#'
#' @param matx A symmetric, positive-definite numeric matrix.
#'
#' @return An upper triangular matrix `R` such that `t(R) %*% R` is equal to
#'   the original matrix `matx`.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' A <- matrix(c(4, 12, -16, 12, 37, -43, -16, -43, 98), nrow = 3)
#' R <- dpdch(A)
#' print(R)
#' #      [,1] [,2] [,3]
#' # [1,]    2    6   -8
#' # [2,]    0    1    5
#' # [3,]    0    0    3
#'
#' # Verify the decomposition: t(R) %*% R should equal A
#' print(t(R) %*% R)
dpdch <- function(matx) {
  # R's chol() function directly computes the Cholesky decomposition.
  # It returns the upper triangular factor R.
  return(chol(matx))
}

#' Solve a Lower Triangular System (Forward Substitution)
#'
#' @description
#' Solves the linear system `L %*% x = b` for `x`, where `L` is a lower
#' triangular matrix. This is an R equivalent of the `subfd` function.
#'
#' @details
#' This function is a wrapper for R's built-in `forwardsolve()` function,
#' which is optimized for this task. The C version modifies the `b` vector
#' in place, while this R function returns the solution vector `x`.
#'
#' @param l A lower triangular numeric matrix.
#' @param b A numeric vector representing the right-hand side of the equation.
#'
#' @return A numeric vector `x` which is the solution to the system.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' L <- matrix(c(2, 0, 0, 6, 1, 0, -8, 5, 3), nrow = 3)
#' b <- c(4, 18, -2)
#' x <- subfd(L, b)
#' print(x) # Expected: c(2, 6, 8)
subfd <- function(l, b) {
  return(forwardsolve(l, b))
}

#' Solve an Upper Triangular System (Backward Substitution)
#'
#' @description
#' Solves the linear system `U %*% x = b` for `x`, where `U` is an upper
#' triangular matrix. This is an R equivalent of the `subbk` function.
#'
#' @details
#' This function is a wrapper for R's built-in `backsolve()` function,
#' which is optimized for this task. The C version modifies the `b` vector
#' in place, while this R function returns the solution vector `x`.
#'
#' @param u An upper triangular numeric matrix.
#' @param b A numeric vector representing the right-hand side of the equation.
#'
#' @return A numeric vector `x` which is the solution to the system.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' U <- matrix(c(2, 6, -8, 0, 1, 5, 0, 0, 3), nrow = 3)
#' b <- c(4, 18, -2) # This is a different 'b' than the subfd example
#' x <- subbk(U, b)
#' print(x)
subbk <- function(u, b) {
  return(backsolve(u, b))
}

#' Solve a Linear System using Cholesky Decomposition
#'
#' @description
#' Solves the linear system `a %*% x = b` for `x`, where `a` is a symmetric,
#' positive-definite matrix. It uses Cholesky decomposition (`a = R'R`),
#' followed by forward and backward substitution. This function combines the
#' logic of `dpdch`, `subfd`, and `subbk`.
#'
#' @details
#' The process involves two steps:
#' 1. Solve `R' * y = b` for `y` using forward substitution.
#' 2. Solve `R * x = y` for `x` using backward substitution.
#'
#' This is a more direct R translation of the `chsol` C function. For general
#' use in R, `solve(a, b)` is often more direct.
#'
#' @param a A symmetric, positive-definite numeric matrix.
#' @param b A numeric vector.
#'
#' @return A numeric vector `x` which is the solution to the system.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' A <- matrix(c(4, 12, -16, 12, 37, -43, -16, -43, 98), nrow = 3)
#' b_vec <- c(4, 1, 2)
#' x_sol <- chsol(A, b_vec)
#' print(x_sol)
#' # Verify the solution: A %*% x_sol should be close to b_vec
#' print(A %*% x_sol)
chsol <- function(a, b) {
  # Step 1: Perform Cholesky decomposition to get R (upper triangular)
  R <- dpdch(a)

  # Step 2: Solve R' * y = b for y using forward substitution
  # t(R) gives the lower triangular matrix L
  y <- subfd(t(R), b)

  # Step 3: Solve R * x = y for x using backward substitution
  x <- subbk(R, y)

  return(x)
}

#' Evaluate a Linear Polynomial (Linear Interpolation)
#'
#' @description
#' Evaluates a linear function passing through two points `(x0, y0)` and
#' `(x1, y1)` at a new point `xvalue`. This is effectively linear interpolation.
#' This function is an R equivalent of `linearPoly` from the C library.
#'
#' @details
#' This function uses R's built-in `stats::approxfun` to create an
#' interpolation function and then applies it to the `xvalue`. This is more
#' robust and idiomatic than a direct translation of the C code's slope-intercept
#' calculation. The function can handle a vector of `xvalue`s.
#'
#' @param x0 A numeric value for the x-coordinate of the first point.
#' @param y0 A numeric value for the y-coordinate of the first point.
#' @param x1 A numeric value for the x-coordinate of the second point.
#' @param y1 A numeric value for the y-coordinate of the second point.
#' @param xvalue A numeric vector of x-values at which to evaluate the function.
#'
#' @return A numeric vector of the interpolated y-values.
#' @author Jaehoon Seol (Original C code), Google's Gemini (R translation)
#'
#' @examples
#' # Find the y-value at x=5 for a line passing through (0,0) and (10,20)
#' linear_poly(x0 = 0, y0 = 0, x1 = 10, y1 = 20, xvalue = 5)
#' #> [1] 10
#'
#' # Can also handle vectors
#' linear_poly(x0 = 0, y0 = 0, x1 = 10, y1 = 20, xvalue = c(2.5, 5, 7.5))
#' #> [1]  5 10 15
linear_poly <- function(x0, y0, x1, y1, xvalue) {
  # Check for a near-zero slope, mimicking the C function's warning
  if (abs(x1 - x0) < .Machine$double.eps) {
    if (abs(y1 - y0) > .Machine$double.eps) {
      warning("Vertical line: returning mean of y0 and y1")
      return(mean(y0, y1))
    }
  }

  # Create a linear interpolation function from the two points
  interp_fun <- stats::approxfun(x = c(x0, x1), y = c(y0, y1), rule = 2)

  # Apply the function to the new x-value(s) and return the result
  return(interp_fun(xvalue))
}

#' Evaluate a Cubic Polynomial
#'
#' @description
#' Evaluates a cubic polynomial of the form:
#' `f(x) = a + b*(x-x_left) + c*(x-x_left)^2 + d*(x-x_left)^3`
#' for a given `xvalue`.
#'
#' @param x_left The left boundary of the interval for the polynomial.
#' @param a The constant coefficient.
#' @param b The linear coefficient.
#' @param c The quadratic coefficient.
#' @param d The cubic coefficient.
#' @param xvalue The x-value at which to evaluate the polynomial.
#'
#' @return The numeric value of the polynomial at `xvalue`.
cubic_poly <- function(x_left, a, b, c, d, xvalue) {
  temp_x <- xvalue - x_left
  return(a + b * temp_x + c * temp_x^2 + d * temp_x^3)
}

#' Calculate Coefficients for a Smoothing Cubic Spline
#'
#' @description
#' This is the core function for calculating the coefficients of a smoothing
#' cubic spline. It implements the algorithm described by Reinsch (1967) and
#' referenced in Kolen & Brennan (2004).
#'
#' @details
#' The function iteratively finds an optimal smoothing parameter `p` and solves
#' a system of linear equations to determine the spline coefficients.
#'
#' @param x A numeric vector of the x-coordinates of the data points.
#' @param y A numeric vector of the y-coordinates of the data points.
#' @param dyi A numeric vector of the standard deviations of the y-coordinates.
#' @param s A fidelity constant that controls the closeness of the spline to the data.
#'
#' @return A matrix containing the spline coefficients (a, b, c, d).
sspline <- function(x, y, dyi, s) {
  n <- length(x) - 1
  s <- s * (x[n + 1] - x[1] + 1)

  h <- diff(x)

  # Setup matrices T and Q
  T_mat <- diag(2 * (h[-n] + h[-1]) / 3)
  diag(T_mat[-1, ]) <- h[-c(1, n)] / 3
  diag(T_mat[, -1]) <- h[-c(1, n)] / 3

  Q <- matrix(0, nrow = n - 1, ncol = n + 1)
  for (i in 1:(n - 1)) {
    Q[i, i] <- 1 / h[i]
    Q[i, i + 1] <- - (1 / h[i] + 1 / h[i + 1])
    Q[i, i + 2] <- 1 / h[i + 1]
  }

  y2 <- Q %*% y

  p <- 0.0
  np <- 1.0

  ii <- 0
  while (abs(np - p) > 1e-9 && ii < 35) {
    p <- np
    # CORRECTED FORMULA: Q %*% D^2 %*% t(Q)
    mpt2 <- Q %*% diag(dyi^2) %*% t(Q) + p * T_mat

    u <- chsol(mpt2, y2)
    v <- diag(dyi) %*% t(Q) %*% u
    e <- sum(v^2)

    Tu <- T_mat %*% u
    w <- chsol(mpt2, Tu)
    Tw <- T_mat %*% w

    f_val <- sum(u * Tu)
    g_val <- sum(u * Tw)

    if (s != 0) {
      if ((f_val - p * g_val) == 0) { # Avoid division by zero
        break
      }
      np <- p + sqrt(e / s) * (e - sqrt(s * e)) / (f_val - p * g_val)
    } else {
      break # s=0 implies interpolation, no need to iterate p
    }

    if (e <= s) break
    ii <- ii + 1
  }

  a <- y - dyi^2 * (t(Q) %*% (p * u))
  c_coef <- c(0, p * u, 0) # Pad with 0 at both ends

  d_coef <- diff(c_coef) / (3 * h)
  b_coef <- (diff(a) / h) - (c_coef[-(n + 1)] + d_coef * h) * h

  return(cbind(a = a[-(n+1)], b = b_coef[,1], c = c_coef[-(n+1)], d = d_coef))
}
#' Apply Post-Smoothing Using Cubic Splines
#'
#' @description
#' Implements the post-smoothing method described in Kolen & Brennan (2004),
#' which involves fitting a cubic spline to a portion of the data and using
#' linear interpolation for the tails.
#'
#' @param xvalues A numeric vector of raw scores.
#' @param yvalues A numeric vector of equated scores corresponding to `xvalues`.
#' @param dyi A numeric vector of standard errors for the `yvalues`.
#' @param s The fidelity (smoothing) constant for the spline.
#' @param xlow The index of the lowest score to be included in the spline fit.
#' @param xhigh The index of the highest score to be included in the spline fit.
#' @param ky The maximum possible score on the target scale (Y).
#' @param vectX The x-coordinates where the smoothed function should be evaluated.
#'
#' @return A list containing the smoothed y-values (`vectY`) and the spline
#'   coefficient matrix (`cmat`).
post_smooth <- function(xvalues, yvalues, dyi, s, xlow, xhigh, ky, vectX) {

  num_spline_pts <- xhigh - xlow + 1

  # R uses 1-based indexing, C uses 0-based. Adjust indices.
  xlow_r <- xlow + 1
  xhigh_r <- xhigh + 1

  cmat <- sspline(
    x = xvalues[xlow_r:xhigh_r],
    y = yvalues[xlow_r:xhigh_r],
    dyi = dyi[xlow_r:xhigh_r],
    s = s
  )

  # Evaluate the spline at the knots to get the smoothed values
  dy_xlow <- cubic_poly(xvalues[xlow_r], cmat[1, "a"], cmat[1, "b"], cmat[1, "c"], cmat[1, "d"], xvalues[xlow_r])
  dy_xhigh <- cubic_poly(xvalues[xhigh_r - 1], cmat[nrow(cmat), "a"], cmat[nrow(cmat), "b"], cmat[nrow(cmat), "c"], cmat[nrow(cmat), "d"], xvalues[xhigh_r])

  kx <- tail(xvalues, 1)

  # Interpolate/evaluate for each point in vectX
  vectY <- sapply(vectX, function(xv) {
    if (xv < xvalues[xlow_r]) {
      linear_poly(-0.5, -0.5, xvalues[xlow_r], dy_xlow, xv)
    } else if (xv > xvalues[xhigh_r]) {
      linear_poly(xvalues[xhigh_r], dy_xhigh, kx + 0.5, ky + 0.5, xv)
    } else {
      j <- findInterval(xv, xvalues[xlow_r:xhigh_r], rightmost.closed = TRUE)
      cubic_poly(xvalues[xlow_r + j - 1], cmat[j, "a"], cmat[j, "b"], cmat[j, "c"], cmat[j, "d"], xv)
    }
  })

  return(list(vectY = vectY, cmat = cmat))
}

#' Find the Inverse of a Cubic Polynomial
#'
#' @description
#' Finds the x-value for a given y-value for a monotonic cubic polynomial,
#' effectively computing the inverse.
#'
#' @details
#' This function uses R's `uniroot` function to find the root of the equation
#' `f(x) - y = 0`, which is a robust and efficient method.
#'
#' @param yvalue The y-value for which to find the corresponding x-value.
#' @param x_left The left boundary of the interval for the polynomial.
#' @param a,b,c,d The coefficients of the cubic polynomial.
#'
#' @return The x-value corresponding to `yvalue`.
inverse_cubic_poly <- function(yvalue, x_left, a, b, c, d) {
  # Define the function whose root we want to find: f(x) - y = 0
  f_to_solve <- function(x) {
    cubic_poly(x_left, a, b, c, d, x) - yvalue
  }

  # Use uniroot to find the solution within the interval
  # The interval is [x_left, x_left + 1] assuming an increment of 1
  # A more robust interval might be needed depending on the spline's behavior
  root_result <- try(uniroot(f_to_solve, interval = c(x_left, x_left + 1)), silent = TRUE)

  if (inherits(root_result, "try-error")) {
    # Handle cases where the root is not in the interval, maybe expand search
    return(NA)
  } else {
    return(root_result$root)
  }
}

#' Compute the Inverse of a Full Cubic Spline Function
#'
#' @description
#' Computes the inverse of a complete cubic spline function at specified points.
#'
#' @param ynodes The original y-nodes (knots) of the spline.
#' @param cmat The coefficient matrix of the spline.
#' @param vectX The x-values at which to evaluate the inverse function.
#'
#' @return A vector of y-values from the inverse function.
inverse_sspline <- function(ynodes, cmat, vectX) {
  # First, find the x-values at the original y-nodes
  xnodes <- sapply(1:nrow(cmat), function(i) {
    cubic_poly(ynodes[i], cmat[i, "a"], cmat[i, "b"], cmat[i, "c"], cmat[i, "d"], ynodes[i])
  })
  xnodes <- c(xnodes, tail(yvalues,1)) # Add the final point

  # Now, for each value in vectX, find the corresponding inverse y-value
  sapply(vectX, function(xv) {
    # Find which segment of the spline xv falls into
    j <- findInterval(xv, xnodes, rightmost.closed = TRUE)
    if (j == 0) j <- 1 # Handle edge case

    inverse_cubic_poly(
      yvalue = xv,
      x_left = ynodes[j],
      a = cmat[j, "a"],
      b = cmat[j, "b"],
      c = cmat[j, "c"],
      d = cmat[j, "d"]
    )
  })
}

#' Compute the Inverse of the Post-Smoothed Equating Function
#'
#' @description
#' Computes the inverse of the complete post-smoothed equating function,
#' including the linear tail sections.
#'
#' @param yvalues Original y-scores (knots).
#' @param xvalues Original equated x-scores.
#' @param dxi Standard errors of equated scores.
#' @param s Smoothing parameter.
#' @param ylow,yhigh Indices for the spline fitting range.
#' @param kx Maximum score on the X scale.
#' @param vectX The x-values at which to evaluate the inverse.
#'
#' @return A vector of y-values from the inverse function.
inverse_post_smooth <- function(yvalues, xvalues, dxi, s, ylow, yhigh, kx, vectX) {

  num_spline_pts <- yhigh - ylow + 1
  ylow_r <- ylow + 1
  yhigh_r <- yhigh + 1

  # Fit the spline d_x(y)
  cmat <- sspline(
    x = yvalues[ylow_r:yhigh_r],
    y = xvalues[ylow_r:yhigh_r],
    dyi = dxi[ylow_r:yhigh_r],
    s = s
  )

  # Evaluate the spline at the knots to get the smoothed values
  dx_ylow <- cubic_poly(yvalues[ylow_r], cmat[1, "a"], cmat[1, "b"], cmat[1, "c"], cmat[1, "d"], yvalues[ylow_r])
  dx_yhigh <- cubic_poly(yvalues[yhigh_r - 1], cmat[nrow(cmat), "a"], cmat[nrow(cmat), "b"], cmat[nrow(cmat), "c"], cmat[nrow(cmat), "d"], yvalues[yhigh_r])

  ky <- tail(yvalues, 1)

  # Evaluate the inverse for each point in vectX
  sapply(vectX, function(xv) {
    if (xv < dx_ylow) {
      linear_poly(x0 = -0.5, y0 = -0.5, x1 = dx_ylow, y1 = yvalues[ylow_r], xvalue = xv)
    } else if (xv > dx_yhigh) {
      linear_poly(x0 = dx_yhigh, y0 = yvalues[yhigh_r], x1 = kx + 0.5, y1 = ky + 0.5, xvalue = xv)
    } else {
      inverse_sspline(ynodes = yvalues[ylow_r:yhigh_r], cmat = cmat, vectX = xv)
    }
  })
}
