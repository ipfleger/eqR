
run_equating <- function(equating_recipe){
  # For each row in plan
    # For each equating method and type and smoothing
      # 1. gather the appropriate information
      # 2. call the appropriate wrapper function
      # 3. Store the results
  # return the results
}


# How about communicating and evaluating the results?







eq <- init_equating() |>
  add_form(ctabs_a, name = "Form A", id_cols = 'id') |>
  add_form(ctabs_b, name = "Form B", id_cols = 'id') |>
  # add_form(ctabs_c, name = "Form C", id_cols = 'id') |>
  add_plan("Form B" = "Form A") |>
  add_design("NE") |>
  add_method(
    method = "equipercentile",
    smooth = "continuized_log_linear",
    type = "chained"
  )
anchors <-  get_anchors(eqb)[[1]][1:5]
k <- 5
fa <- sample(colnames(ctabs_a)[!colnames(ctabs_a) %in% anchors][-1], k)
fb <- sample(colnames(ctabs_b)[!colnames(ctabs_b) %in% anchors][-1], k)


eq <- init_equating() |>
  add_form(ctabs_a[,c(fa, anchors)], name = "Form A", id_cols = 'id') |>
  add_form(ctabs_b[,c(fb, anchors)], name = "Form B", id_cols = 'id') |>
  # add_form(ctabs_c, name = "Form C", id_cols = 'id') |>
  add_plan("Form B" = "Form A") |>
  add_design("NE") |>
  add_method(
    method = "equipercentile",
    smooth = "continuized_log_linear",
    type = "chained"
  )


blls <- lapply(str_split(rownames(eq@plan[1]), pattern = ";")[[1]] |> `names<-`(str_split(rownames(eq@plan[1]), pattern = ";")[[1]]), \(frm) {

  points <- apply(eq@data[[frm]][colnames(eq@data[[frm]])], 2, max, na.rm = TRUE)
  u_points <- sum(points[!colnames(eq@data[[frm]]) %in% get_anchors(eq)[[rownames(eq@plan[1])]]])
  v_points <- sum(points[colnames(eq@data[[frm]]) %in% get_anchors(eq)[[rownames(eq@plan[1])]]])
  u_score = rowSums(eq@data[[frm]][!colnames(eq@data[[frm]]) %in% get_anchors(eq)[[rownames(eq@plan[1])]]], na.rm = TRUE)

  v_score = rowSums(eq@data[[frm]][colnames(eq@data[[frm]]) %in% get_anchors(eq)[[rownames(eq@plan[1])]]], na.rm = TRUE)
  rv = cor(u_score, v_score, use = "p")
  freq_tab <- table(u_score = factor(u_score, levels = 0:u_points), v_score = factor(v_score, levels = 0:v_points))

  so <- eq@methods$`CG E C Z`$options
  n <- nrow(eq@data[[frm]])
  nsu <- u_points + 1
  min_u <- so$min_u
  incu <- so$inc_u
  nsv <- v_points + 1
  minv <- so$min_v
  incv <-so$inc_v
  nct <- as.vector(t(freq_tab))
  scale = eq@methods$`CG E C Z`$smooth != "Z" # Set to FALSE for continuized_log_linear, TRUE for all others


  # There should be something here to verify that the smoothing fits the model well.

  smooth_bll(n = n, nsu = nsu, minu = min_u, incu = incu, nsv = nsv,
             minv = minv, incv = incv, nct = nct, cu = so$cu, cv = so$cv,
             cuv = so$cuv, cpm = so$cpm, scale = scale, crit = so$crit)

})

  equated_scores <- equate_cll(
    design = "NEAT_CHN",
    bivar1 = blls[[1]],
    bivar2 = blls[[2]]
  )

  R = 10

  system.time({
    boot_cll <- lapply(1:R,\(i) {

      bv1 <- smooth_bll(n = blls[[1]]$n,
                        nsu = blls[[1]]$nsx,
                        minu = blls[[1]]$minx,
                        incu = blls[[1]]$incx,
                        nsv = blls[[1]]$nsv,
                        minv = blls[[1]]$minv,
                        incv = blls[[1]]$incv,
                        nct = as.vector(t(bootstrap_parametric_biv(n = blls[[1]]$n, smoothed_bfd = matrix(blls[[1]]$mct/blls[[1]]$n, nrow = blls[[1]]$nsx, ncol = blls[[1]]$nsv)))),
                        cu = blls[[1]]$cu,
                        cv = blls[[1]]$cv,
                        cuv = blls[[1]]$cuv,
                        cpm = blls[[1]]$cpm,
                        scale = eq@methods$`CG E C Z`$smooth != "Z",
                        crit = blls[[1]]$crit)

      bv2 <- smooth_bll(n = blls[[2]]$n,
                        nsu = blls[[2]]$nsx,
                        minu = blls[[2]]$minx,
                        incu = blls[[2]]$incx,
                        nsv = blls[[2]]$nsv,
                        minv = blls[[2]]$minv,
                        incv = blls[[2]]$incv,
                        nct = as.vector(t(bootstrap_parametric_biv(n = blls[[2]]$n, smoothed_bfd = matrix(blls[[2]]$mct/blls[[2]]$n, nrow = blls[[2]]$nsx, ncol = blls[[2]]$nsv)))),
                        cu = blls[[2]]$cu,
                        cv = blls[[2]]$cv,
                        cuv = blls[[2]]$cuv,
                        cpm = blls[[2]]$cpm,
                        scale = eq@methods$`CG E C Z`$smooth != "Z",
                        crit = blls[[2]]$crit)

      equate_cll(
        design = "NEAT_CHN",
        bivar1 = bv1,
        bivar2 = bv2,
      )

    })

  })




boot_results <- do.call(cbind, boot_cll)


# ggplot2 is easier to style, relative is faster.


 round(equated_scores, 2)

# 1. Model and Design Specification





 print_cll_cg(title = "test",
              pdata = blls,
              conversion_table = conversion_table)

basetheme::basetheme("minimal")
plot_equivalent(conversion_table, relative = TRUE)
# The least amount of ink required to communi


# 2. Equated Raw Scores

 conversion_table <- data.frame(score = 1:blls[[2]]$nsx-1,
                                equivalent = rowMeans(boot_results),
                                equivalent_no_bs = equated_scores,
                                se = apply(boot_results, 1, sd),
                                lower_bound_95 = apply(boot_results, 1, quantile, .025),
                                lower_bound_50 = apply(boot_results, 1, quantile, .25),
                                upper_bound_50 = apply(boot_results, 1, quantile, .75),
                                upper_bound_95 = apply(boot_results, 1, quantile, .975))



# 3. Moments of the Equated Scores


  # Correct way to create the long data frame
  long_df <- data.frame(
    u_score = rep(0:(blls[[1]]$nsx - 1), times = blls[[1]]$nsv),
    v_score = rep(0:(blls[[1]]$nsv - 1), each = blls[[1]]$nsx),
    frequency = blls[[1]]$mct
  )
  fig <- plotly::plot_ly(
    data = long_df,
    x = ~u_score,
    y = ~v_score,
    z = ~frequency,
    type = 'mesh3d',
    intensity = ~frequency,
    colorscale = 'Viridis'
  )
  fig <- fig %>% plotly::layout(
    title = "Smoothed Bivariate Score Distribution",
    scene = list(
      xaxis = list(title = "Score on Unique Items (u)"),
      yaxis = list(title = "Score on Anchor Items (v)"),
      zaxis = list(title = "Smoothed Frequency")
    )
  )
  print(fig)

# Garbage ---------------------
  args <- list(...) # or is args defaults passed through already, like increment, minimum, min_anch, inc_anch,
  # For each form
  lapply(str_split(rownames(eq@plan[1]), pattern = ";")[[1]], \(frm) {
    lapply(names(eq@methods), \(title){
      categories <- str_split(title, pattern = " ")[[1]]
      equate(design = categories[1], method = categories[3], smooth = categories[4], args)
    })
  })


# if conditions for prepare_equate_cll are met,
  # bivar_list <- prepare_equate_cll
  # equate_cll(bivar_list)

 # Sensible defaults -----------------------------------------------------------
# under_add_methods we may need to break out and have wrappers to the wrappers.
# not so much that as a function to prepare the data for the given wrapper.
# Inside add_method, if smooth is log-linear and no degrees are given...



# Information for the add_method details ---------------------------------------
# I'm not thrilled with this layout. I'd prefer something more indepth.
# I'd also like a section sharing what the equivalents would be to the original code.
# This would probably be under examples, I'd have a Wrapper_CC with the equivalent add_method.
# Method
# linear (and mean)
## type
## smoothing: there are no smoothing methods for linear equating
# equipercentile
## type
## smoothing

#' @param equate_recipe An object of class `equate_recipe`.
#' @param method A character string specifying the primary equating method.
#' @param smooth A character string specifying the smoothing procedure. Defaults to `"none"`.
#' @param type A character string specifying the calculation type or sub-method. Defaults to the most common type for the selected method.
#' @param ... Additional named arguments passed to the smoothing function.
#'
#' @details
#' This function adds a specific equating method to be performed. The required
#' arguments under `...` depend on the combination of `method` and `smooth` chosen.
#'
#' ### Method: `equipercentile`
#'
#' Equipercentile equating is a process that aligns the score distributions of
#' two test forms. Unlike linear equating, which only considers the mean and
#' standard deviation, equipercentile equating matches the entire score
#' distribution by finding scores on the two forms that have the same
#' percentile rank. This often results in a more accurate, non-linear
#' relationship between the score scales.
#'
#' ### Smoothing: `continuized_log_linear`
#'
#' This method first fits a log-linear model to the discrete score distribution(s)
#' and then uses the resulting model parameters to define a smooth, continuous
#' probability density function. This continuous function is then used in the
#' equipercentile equating process. It is controlled by the following parameters,
#' which are passed via `...`:
#'
#' * `degree_u`: The polynomial degree for the first variable (e.g., Form X).
#'     Defaults to `6`.
#' * `degree_v`: The polynomial degree for the second variable (e.g., the anchor
#'     test V). Defaults to `6`.
#' * `degree_uv`: The number of cross-product terms to model the interaction.
#'     Defaults to `1`.
#' * `cpm`: A matrix specifying the powers for the cross-product terms.
#'     Defaults to `matrix(c(1, 1), nrow = 1)`, which corresponds to the
#'     term $u^1v^1$.
#'
#' These parameters define the formula used for the log-linear model, which you
#' can see in the "Log-Linear Model Formula" document.
#'
#' ### Type: `chained`
#'
#' This `type` is used with the Common-Item Non-Equivalent Groups (`CG`) design.
#' It performs the equating in two steps:
#' 1.  It equates the "from" form (e.g., Form B) to the scale of the anchor test (V)
#'     using the data from the group that took Form B.
#' 2.  It then equates those anchor-scale scores to the scale of the "to" form
#'     (e.g., Form A) using the data from the group that took Form A.
#'
#' This method does not require creating a synthetic population.
#'
#' ### Calculation Types (`type`) for the CG Design
#'
#' The `type` argument specifies a sub-method, primarily for the Common-Item (`CG`)
#' design when `method = "equipercentile"`. The underlying C library wrappers
#' (e.g., `Wrapper_CN`, `Wrapper_CL`) use single-letter codes to specify the
#' calculation type [cite: p. 69]. This R function maps the following
#' user-friendly names to those codes:
#'
#' * `"frequency"` (default): Runs Frequency Estimation (FE) and Braun-Holland under FE. (Code: `E`)
#' * `"chained"`: Runs Chained Equipercentile equating. (Code: `C`)
#' * `"modified_frequency"`: Runs Modified Frequency Estimation (MFE) and Braun-Holland under MFE. (Code: `F`)
#' * `"all"`: Runs all available methods (FE, MFE, Chained, and associated Braun-Holland). (Code: `A`)
#' * `"G"`: Runs FE, BH-FE, MFE, and BH-MFE.
#' * `"H"`: Runs FE, BH-FE, and Chained.
#'
