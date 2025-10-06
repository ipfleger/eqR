#' Standardize IRT Item Parameters
#'
#' This function takes item parameters from various sources and converts them
#' into a standardized data.frame format with slope, intercept, and guess columns.
#'
#' @param irt_pars An object containing IRT parameters. Can be a data.frame,
#'   a mirt object, a TAM object, or a list of these.
#' @return A list of standardized item parameter data.frames.
#' @keywords internal
irt_coefs <- function(irt_pars){
  if(is.list(irt_pars) && !is.data.frame(irt_pars)){
    return(lapply(irt_pars, irt_coefs))
  }

  if(is.data.frame(irt_pars)){
    # Standardize from a generic data frame
    # Handles both a/b/c and slope/intercept/guess formats
    df <- data.frame(
      item = irt_pars$item,
      slope = irt_pars$slope %||% irt_pars$a,
      intercept = irt_pars$intercept %||% (-1 * (irt_pars$a %||% 0) * (irt_pars$b %||% 0)),
      guess = irt_pars$guess %||% irt_pars$c %||% 0
    )
    return(df)

  } else if((attr(class(irt_pars), "package") %||% "not") == "mirt"){
    # Extract from a mirt model object
    if (!requireNamespace("mirt", quietly = TRUE)) {
      cli::cli_abort("The 'mirt' package is required to extract coefficients from a mirt object.")
    }
    coefs <- mirt::coef(irt_pars, simplify = TRUE)$items
    df <- data.frame(
      item = rownames(coefs),
      slope = coefs[, "a1"],
      intercept = coefs[, "d"],
      guess = coefs[, "g"]
    )
    return(df)

  } else if (any(grepl("tam", class(irt_pars)))) {
    # Extract from a TAM model object
    df <- irt_pars$item[, c("item", "B.Cat1.Dim1", "AXsi_.Cat1", "guess")]
    colnames(df) <- c("item", "slope", "intercept", "guess")
    return(df)

  } else {
    cli::cli_abort("The format of 'item_params' is not recognized. Please provide a data.frame, mirt object, TAM object, or a list of these.")
  }
}
