# Single Group
pse <- readRDS("~/statwise ecosystem 2.0/lyzer/data/pse_google.rds")

forms <- formr:::build(is = bl@is |> dplyr::mutate(decision = "Use"),
                       f = 3,
                       ipo = bl@ipo[,c("obj", "per")] |> dplyr::mutate(per = round(per/2.86)),
                       okol = data.frame(forms = c("form1,form2"), ol = .2))

bl <- formr:::pse2bl(pse)
bl@forms <- forms |> dplyr::rename(`Form A` = form1,
                            `Form B` = form2,
                            `Form C` = form3)

devtools::load_all(".")
ctabs_a <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form A`), which = "s", id_cols = c("id"))
ctabs_b <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form B`), which = "s", id_cols = c("id"))
ctabs_c <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form C`), which = "s", id_cols = c("id"))

eq <- init_equating() |>
  add_form(ctabs_a, name = "Form A", id_cols = 'id', min_score = 0, max_score = 39) |>
  add_form(ctabs_b, name = "Form B", id_cols = 'id', min_score = 0, max_score = 39) |>
  add_form(ctabs_c, name = "Form C", id_cols = 'id', min_score = 0, max_score = 39) |>
  add_plan(`Form A` ~ `Form C` + `Form B`)


# Single Group Design -----------------------------------
# Mean equating
single_mean <- eq |>
  add_design("SG") |> # Add the design
  add_method(
    method = "linear", mean_only = TRUE # Mean equating
  ) |>
  run_equating(boot_replications = 1000, boot_type = "norm")
# Linear equating
single_linear <- eq |>
  add_design("SG") |> # Add the design
  add_method(
    method = "linear",
  ) |>
  run_equating(boot_replications = 1000, boot_type = "basic")
# Equipercentile equating
single_equipercentile <- eq |>
  add_design("SG") |>
  add_method(
    method = "equipercentile" # equipercentile equating no smoothing
  ) |>
  run_equating(boot_replications = 1000) # bootstrap options are from boot::boot

single_equip_loglinear <- eq |>
  add_design("SG") |>
  add_method(
    method = "equipercentile", smooth = "log_linear" # equipercentile equating no smoothing
  ) |>
  run_equating(boot_replications = 1000)

single_equip_continuized_loglinear <- eq |>
  add_design("SG") |>
  add_method(
    method = "equipercentile", smooth = "continuized_log_linear" # equipercentile equating no smoothing
  ) |>
  run_equating(boot_replications = 1)  # Don't run this with 1000 boot_replications. It will run until you die. Really any replications will make it take a long time.


plot_equivalent(results = single_equip_continuized_loglinear$`Form C;Form A`$`S E I Z`$`Equipercentile`)
# Missing confidence intervals don't plot well.


# A test:
dat <- data.frame(do.call(cbind, lapply(forms, \(frm) rowSums(eq@data[[frm]], na.rm = TRUE))))
names(dat) <- c("x", "y")
data = dat; i = 1:nrow(dat)
n = equate_none_stat_fun(data, i, score_params, method_options) # is this mfe or fe or bh?
b = equate_bb_stat_fun(data, i, score_params, method_options)
l = equate_loglinear_stat_fun(data, i, score_params, method_options)
# s = equate_spline_stat_fun(data, i, score_params, method_options) # doesn't run
k = equate_kernel_stat_fun(data, i, score_params, method_options) # doesn't seem to be working correctly.
z = equate_cll_sg_stat_fun(data, i, score_params, method_options) # Possibly needs our adjustment
data.frame(x = 0:70, n = n, b = b, l = l, #s = s,
           k = k, z = z) |>
  tidyr::pivot_longer(-x) |>
  ggplot2::ggplot(ggplot2::aes(x = x, y = value, color = name)) +
  ggplot2::geom_line()

# Random Group Design -----------------------------------
# equpercentile equating with smoothing
single_equip_bb <- eq |> add_design("RG") |>
  add_method(
    method = "equipercentile", smooth = "beta_binomial"
  ) |>
  run_equating(boot_replications = 1000)



single$`Form C;Form A`$`S L mean N mean_only`$Mean
plot_equivalent(results = single$`Form C;Form A`$`S E identity N`$`Equipercentile (No Smoothing)`)
plot_equivalent(results = single$`Form B;Form A`$`S E identity N`$`Equipercentile (No Smoothing)`)
single$`Form B;Form A`$`S L linear N`$Linear$plots$score_conversion(relative = FALSE) # I actually really like this though, this feels nice. Darn it, I might actually like the pythonic way of doing it better.
single$`Form C;Form A`$`S L linear N`$Linear$plots$score_conversion(relative = TRUE)

# Common-Item Group Design -----------------------------------


common <- eq |>
  add_design("CG") |>
  # add_method(
  #   method = "linear", type = "all", internal_anchors = TRUE, mean_only = FALSE, w1 = 1 # Can't be zero, synthetic variance is non-positive. What is that about?
  # )|>
  add_method(
    method = "equipercentile", type = c("E"), #smooth = "continuized_log_linear"
  ) |>
  run_equating(boot_type = "bca") # this breaks as a bca | 1: In norm.inter(t, adj.alpha) : extreme order statistics used as endpoints

plot_equivalent(results = common$`Form C;Form A`$`CG E E N`$FrequencyEstimation)
plot_equivalent(results = common$`Form C;Form A`$`CG E E N`$BraunHolland_FE)

# FE
# BH
# MFE
# bca not working for equipercentile


plot_equivalent(results = common$`Form C;Form A`$`CG E C N`$Chained)

common$`Form C;Form A`$`CG E H Z`$FrequencyEstimation
plot_equivalent(results = common$`Form C;Form A`$`CG L all N`$`Levine Observed`)
plot_equivalent(results = common$`Form C;Form A`$`CG E H N`$BraunHolland_FE)


plot_equivalent(results = common$`Form C;Form A`$`CG E H Z`$FrequencyEstimation)
plot_equivalent(results = common$`Form C;Form A`$`CG E H Z`$BraunHolland_FE)
plot_equivalent(results = common$`Form C;Form A`$`CG E H Z`$Chained)
plot_equivalent(results = common2$`Form C;Form A`[[1]]$Chained)

# I want to make sure that the smoothing is being applied at all.

common$`Form C;Form A`$`CG L all N`$Tucker$plots$score_conversion() # let's just call this plot and just have the one.
sum <- summarize_equating_results(common)
sum$`Form C;Form A`$conversion_table

sum$`Form C;Form A`$parameter_comparison
sum1$`Form C;Form A`$parameter_comparison
sum$`Form C;Form A`$moment_comparison

sum$`Form B;Form A`$parameter_comparison
sum$`Form C;Form A`$conversion_plot_gg

eq$`Form C;Form A`$`CG L tucker N`$Tucker$plots$score_conversion(relative = TRUE)
eqb$`Form C;Form A`$`CG L tucker N`$Tucker$plots$score_conversion(relative = TRUE)

class(eq)

# I'm not sure that the plotting functions are working from summary.
# Either way, I'd like to move on to the others. This is the general idea of
# final output. I can make it better later.
# Let's move on to updating the functions.

# We should add a note saying that smoothing options are ignored when they aren't applicable.
# I want to add max, minimum, and increment as options to the add_data statement
# I think these are mostly relevant when the data are tabular or when not all possible scores are
# represented.
# The default should be smarter though. The max should be the sum of points, which is the colSums of the max of the points scored per item (i.e. eq@data[[frm]])
# I am going through the equipercentile method to try and figure out if it works and if it doesn't, why doesn't it. I may have put too much faith in the AI for this one, hopefully it doesn't take just as long to fix as it would have to write it. It seems to only want to run braun holland and frequency estimation, it doesn't seem to want to run continuized log linear. This shouldn't be too hard to get it all going, I am fine breaking it down into smaller bits.



eq <- init_equating() |>
  add_form(ctabs_a, name = "Form A", id_cols = 'id') |>
  add_form(ctabs_b, name = "Form B", id_cols = 'id') |>
  add_plan(`Form A` ~ `Form B`) |>
  add_design("CG") |>
  add_method(
    method = "equipercentile", smooth = "continuized_log_linear"#type = "frequency",
  ) |>
  run_equating(boot_type = "perc")


equipercentile <- function(...){cli::cli_abort("equipercentile equating not yet implemented")}
irt <- function(...){cli::cli_abort("irt equating not yet implemented")}
