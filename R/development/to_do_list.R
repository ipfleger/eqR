# Single Group
bl <- readr::read_rds("data/bl_service_now.rds")

ctabs_a <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form A`), which = "s", id_cols = c("id"))
ctabs_b <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form B`), which = "s", id_cols = c("id"))
ctabs_c <- formr::crosstabs(bl@pse |> dplyr::filter(item %in% bl@forms$`Form C`), which = "s", id_cols = c("id"))

eq <- init_equating() |>
  add_form(ctabs_a, name = "Form A", id_cols = 'id', min_score = 0, max_score = 70) |>
  add_form(ctabs_b, name = "Form B", id_cols = 'id', min_score = 0, max_score = 70) |>
  add_form(ctabs_c, name = "Form C", id_cols = 'id', min_score = 0, max_score = 70) |>
  add_plan(`Form A` ~ `Form C` + `Form B`)
# Single Group Design -----------------------------------

 single <- eq |>
  add_design("single") |> # Add the design
  add_method(
    method = "linear", mean_only = TRUE # Mean equating
  ) |>
  add_method(
    method = "linear" # linear equating
  )  |>
  add_method(
    method = "equipercentile" # equipercentile equating no smoothing
  ) |>
  run_equating(boot_replications = 1000, boot_type = "perc") # bootstrap options are from boot::boot


single$`Form C;Form A`$`S L mean N mean_only`$Mean
plot_equivalent(results = single$`Form C;Form A`$`S E identity N`$`Equipercentile (No Smoothing)`)
plot_equivalent(results = single$`Form B;Form A`$`S E identity N`$`Equipercentile (No Smoothing)`)
single$`Form B;Form A`$`S L linear N`$Linear$plots$score_conversion(relative = FALSE) # I actually really like this though, this feels nice. Darn it, I might actually like the pythonic way of doing it better.
single$`Form C;Form A`$`S L linear N`$Linear$plots$score_conversion(relative = TRUE)

# Common-Item Group Design -----------------------------------


common <- eq |>
  add_design("CG") |>
  add_method(
    method = "linear", type = "all", internal_anchors = TRUE, mean_only = FALSE, w1 = 1 # Can't be zero, synthetic variance is non-positive. What is that about?
  )|>
  add_method(
    method = "equipercentile", type = c("H"), #smooth = "continuized_log_linear"
  ) |>
  run_equating(boot_type = "perc") # this breaks as a bca | 1: In norm.inter(t, adj.alpha) : extreme order statistics used as endpoints

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
