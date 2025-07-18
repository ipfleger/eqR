
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
