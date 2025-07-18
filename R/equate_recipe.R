# Suggested S7 object definition
equate_recipe <- S7::new_class(
  "equate_recipe",
  properties = list(
    forms    = S7::class_list,
    data     = S7::class_list,
    data_ids = S7::class_list,
    plan     = S7::class_data.frame,
    design   = S7::class_character,
    methods  = S7::class_list
  )
)
# Constructor to initialize the recipe
init_equating <- function() {
  equate_recipe()
}


#' Add an exam form's item response data to an equate recipe
#'
#' Adds a new exam form to an `equate_recipe` object by storing its item response
#' data, form name, and optional candidate ID columns. It also attaches the
#' score scale properties (min, max, increment) as attributes to the data.
#'
#' @param equate_recipe An object of class `equate_recipe` to which the form should be added.
#' @param crosstabs A data frame containing item-level response data for the form.
#' @param name A character string giving the name of the exam form (e.g., "Form A").
#' @param id_cols Optional. A character vector specifying which columns in `crosstabs`
#' contain candidate identifiers.
#' @param min_score,max_score,inc Optional numeric values to manually specify the
#'   score scale. If `NULL`, they are calculated from the data.
#'
#' @return The updated `equate_recipe` object containing the new form.
#' @export
add_form <- function(equate_recipe, crosstabs, name, id_cols = NULL, min_score = NULL, max_score = NULL, inc = NULL) {

  data_cols <- !colnames(crosstabs) %in% id_cols
  form_data <- crosstabs[, data_cols]

  # Get the score scale attributes
  range_attrs <- get_range(form_data, form_name = name, min_score = min_score, max_score = max_score, inc = inc)

  # Safely attach the attributes to the data frame
  attr(form_data, "min") <- range_attrs$min
  attr(form_data, "max") <- range_attrs$max
  attr(form_data, "inc") <- range_attrs$inc
  attr(form_data, "range") <- range_attrs$range
  attr(form_data, "points") <- range_attrs$points

  equate_recipe@forms[[name]] <- colnames(form_data)
  equate_recipe@data[[name]] <- form_data
  equate_recipe@data_ids[[name]] <- crosstabs[, !data_cols]

  equate_recipe
}


#' Add or replace an equating plan
#'
#' Adds an equating plan to an `equate_recipe` object using an intuitive
#' formula-based syntax. If no plan is supplied, it auto-generates all-to-all
#' mappings.
#'
#' @param equate_recipe An object of class `equate_recipe` containing equating data and form names.
#' @param ... A series of formulas specifying the equating plan.
#'   - The **left-hand side (LHS)** of the formula specifies the reference form
#'     (the "to" form, or Form Y).
#'   - The **right-hand side (RHS)** specifies the forms to be equated to the
#'     reference form (the "from" forms, or Form X). Multiple "from" forms can
#'     be specified using `+`.
#'
#' @return An updated `equate_recipe` object with the new or modified plan.
#' @export
#'
#' @examples
#' recipe <- init_equating() |>
#'   add_form(form_a_data, name = "Form A") |>
#'   add_form(form_b_data, name = "Form B") |>
#'   add_form(form_c_data, name = "Form C")
#'
#' # Example 1: Equate Form B and Form C to the reference Form A.
#' recipe <- add_plan(recipe,
#'   `Form A` ~ `Form B` + `Form C`
#' )
#'
#' # Example 2: Equate Form B to A, and Form C to B.
#' recipe <- add_plan(recipe,
#'   `Form A` ~ `Form B`,
#'   `Form B` ~ `Form C`
#' )
#'
#' # Example 3: Auto-generate all possible pairings.
#' recipe <- add_plan(recipe)
add_plan <- function(equate_recipe, ...) {

  plan_formulas <- list(...)

  if (length(equate_recipe@plan) != 0) {
    cli::cli_alert_warning("Overwriting existing plan.")
  }

  if (length(plan_formulas) == 0) {
    # Auto-generate plan from available form names
    form_names <- names(equate_recipe@forms)

    if (length(form_names) > 1) {
      plan_df <- expand.grid(from = form_names, to = form_names, stringsAsFactors = FALSE)
      plan_df <- subset(plan_df, from != to)
      equate_recipe@plan <- plan_df
    } else {
      cli::cli_alert_warning("Not enough forms to generate plan automatically.")
    }

  } else {
    if (!all(sapply(plan_formulas, inherits, "formula"))) {
      cli::cli_abort("All arguments provided to `...` must be formulas.")
    }

    # Parse the list of formulas into a data frame
    plan_list <- lapply(plan_formulas, function(f) {
      to_form <- all.vars(stats::update(f, . ~ 0))
      from_forms <- all.vars(stats::update(f, 0 ~ .))

      if (length(to_form) != 1) {
        cli::cli_abort("Each formula must have exactly one reference form on the left-hand side.")
      }

      data.frame(
        from = from_forms,
        to = to_form,
        stringsAsFactors = FALSE
      )
    })

    plan_df <- do.call(rbind, plan_list)
    rownames(plan_df) <- apply(plan_df, 1, paste, collapse = ";")
    equate_recipe@plan <- plan_df
  }

  equate_recipe
}




#' Add a study design to an equate recipe
#'
#' Adds a normalized design code ("S", "R", or "CG") to an `equate_recipe` object
#' based on a user-friendly description of the equating design.
#'
#' @param equate_recipe An object of class `equate_recipe`.
#' @param design A character string describing the study design. Accepts common
#' names and abbreviations for single-group ("S"), random-groups ("R"), or
#' common-item/nonequivalent groups ("CG") designs.
#'
#' @return The updated `equate_recipe` object with the normalized `@design` field.
#' @export
#'
#' @examples
#' equate_recipe |> add_design("single")
#' equate_recipe |> add_design("random-groups")
#' equate_recipe |> add_design("common-item nonequivalent")
add_design <- function(equate_recipe, design) {
  design <- tolower(design)

  # Match normalized design codes
  if (design %in% c("single", "single-group","single group", "s", "sg")) {
    equate_recipe@design <- "S"
  } else if (design %in% c("random", "random-groups", "r", "rg")) {
    equate_recipe@design <- "R"
  } else if (design %in% c(
    "common-item", "common", "common-item nonequivalent", "common-item non-equivalent",
    "non-equivalent", "nonequivalent", "ne", "cg", "c"
  )) {
    equate_recipe@design <- "CG"
  } else {
    cli::cli_abort("Unknown design type: '{design}'")
  }

  equate_recipe
}

#' Get Method Options with Defaults
#'
#' Merges user-specified options for an equating method with a master list of
#' default values. User options take precedence.
#'
#' @param method_spec A single method specification list from the `equate_recipe`
#'   object (e.g., `eq@methods[[title]]`).
#'
#' @return A list containing the complete set of options for the method.
#' @noRd
get_method_options <- function(method_spec) {
  # Master list of all possible default options
  all_defaults <- list(
    # General
    mean_only = FALSE,

    # Linear CG
    w1 = 1,
    internal_anchors = TRUE,

    # Cubic Spline
    s = 0.2,
    prlow = 0.5,
    prhigh = 99.5,

    # Beta-Binomial
    nparm = 4,
    rel = NULL, # This was 0.85, but I don't think it should have a default. By default it is calculated from the data as a correlation, it can be overwritten if necessary. I also kind of doubt that this is the way to do it, it may need to be specified as a list or some other means. We may need to experiment, but ultimately I don't like this.

    # Log-Linear (univariate)
    degree = 3,

    # Continuized Log-Linear
    cu = 6,
    cv = 6,
    cuv = 1,
    cpm = matrix(c(1, 1), nrow = 1),
    crit = 1e-5,
    inc_u = 1,
    min_u = 0,
    min_v = 0,
    inc_v = 1
  )

  # Get user-provided options from the recipe
  user_options <- method_spec$options

  # Merge them - user_options will overwrite defaults
  final_options <- all_defaults
  for (name in names(user_options)) {
    final_options[[name]] <- user_options[[name]]
  }

  return(final_options)
}

equate <- function(forms, method, design, type, eq, title,
                   boot_type = "perc", boot_replications = 1000){
  if(method == "L"){
    linear(forms = forms, design = design, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if(method == "E"){
    equipercentile(forms = forms, design = design, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  } else if(method == "IRT"){
    irt(forms = forms, design = design, type = type, eq = eq, title = title, boot_type = boot_type, boot_replications = boot_replications)
  }
}

#' Add an equating method to a recipe
#'
#' This function adds a specific equating method, including its calculation type
#' and smoothing options, to an `equate_recipe` object. You can add multiple
#' methods to a single recipe, and they will all be executed by `run_equating()`.
#'
#' @param equate_recipe An object of class `equate_recipe`.
#' @param method A character string specifying the primary equating method.
#'   Can be `"linear"`, `"equipercentile"`, or `"irt"`.
#' @param smooth A character string specifying the smoothing procedure.
#'   Defaults to `"none"`. See the "Smoothing & Optional Arguments" section for details.
#' @param type A character string specifying the calculation type or sub-method,
#'   which is particularly relevant for the common-item nonequivalent groups
#'   design. Defaults to the most common type for the selected method. See the
#'   "Method and Type Details" section for all available options.
#' @param ... Additional named arguments passed to the method or smoothing function.
#'   See the "Smoothing & Optional Arguments" section for details.
#'
#' @return The updated `equate_recipe` object with the new method added.
#' @export
#'
#' @section Method and Type Details:
#' The `type` argument selects a specific statistical procedure. Its valid values
#' depend on the `method` and the `design` set in the recipe.
#'
#' \strong{`method = "linear"`}
#' \itemize{
#'   \item For \strong{Single Group (`SG`) and Random Groups (`RG`)} designs:
#'   \itemize{
#'     \item `type = "linear"` (Default): Matches the means and standard deviations of the two forms.
#'     \item `type = "mean"`: Matches only the means, equivalent to setting `mean_only = TRUE`.
#'   }
#'   \item For \strong{Common-Item (`CG`)} design:
#'   \itemize{
#'     \item `type = "tucker"` (Default): The classic Tucker method for CINEG.
#'     \item `type = "levine_observed"`: The Levine observed score method.
#'     \item `type = "levine_true"`: The Levine true score method, which requires an internal anchor.
#'     \item `type = "chained"`: Chained linear equating through the anchor test.
#'     \item You can also provide a vector of types, e.g., `type = c("tucker", "chained")`, to run multiple methods at once.
#'   }
#' }
#'
#' \strong{`method = "equipercentile"`}
#' \itemize{
#'   \item For \strong{Single Group (`SG`) and Random Groups (`RG`)} designs, the `type` argument is ignored. The standard equipercentile equating procedure is used.
#'   \item For \strong{Common-Item (`CG`)} design:
#'   \itemize{
#'     \item `type = "frequency"` (Default): The Frequency Estimation method.
#'     \item `type = "chained"`: The Chained Equipercentile method.
#'     \item `type = "modified_frequency"`: The Modified Frequency Estimation (MFE) method, which adjusts for anchor reliability. Also runs Braun-Holland under MFE.
#'   }
#' }
#'
#' @section Smoothing & Optional Arguments:
#' The `smooth` argument selects a procedure to smooth the score distributions.
#' Optional parameters for both smoothing and equating methods can be passed via `...`.
#'
#' \strong{Smoothing Procedures}
#' \itemize{
#'   \item \code{smooth = "none"}: (Default) No smoothing is applied.
#'   \item \code{smooth = "beta_binomial"}: Four-parameter beta-binomial smoothing.
#'     \itemize{
#'       \item \code{nparm}: The number of parameters for the beta distribution (2 or 4). Integer. *Default: \code{4}*.
#'       \item \code{rel}: The reliability of the test (e.g., KR-20). Numeric. *Default: \code{0.85}*.
#'     }
#'   \item \code{smooth = "log_linear"}: Univariate log-linear presmoothing.
#'     \itemize{
#'       \item \code{degree}: The degree of the polynomial to fit. Integer. *Default: \code{3}*.
#'     }
#'   \item \code{smooth = "cubic_spline"}: Post-smoothing using cubic splines.
#'     \itemize{
#'       \item \code{s}: The smoothing parameter. Numeric. *Default: \code{0.2}*.
#'     }
#'   \item \code{smooth = "kernel"}: Kernel equating, which uses a continuized score distribution.
#' }
#'
#' \strong{Linear Method Options}
#' \itemize{
#'   \item \code{mean_only}: Sets the slope to 1 to perform mean equating. Logical. *Default: \code{FALSE}*.
#'   \item \code{w1}: For linear CG methods, the weight given to population 1 in the synthetic population. Numeric. *Default: \code{1}*.
#'   \item \code{internal_anchors}: For linear CG methods, specifies if the anchor test is internal. Logical. *Default: \code{TRUE}*.
#' }
#'
#' @examples
#' recipe <- init_equating() |> add_design("common-item")
#' # Add a linear method that runs both Tucker and Chained equating
#' recipe <- recipe |> add_method("linear", type = c("tucker", "chained"))
#' # Add a mean equating method (slope fixed to 1)
#' recipe <- recipe |> add_method("linear", type = "tucker", mean_only = TRUE)
add_method <- function(equate_recipe, method, type = "default", smooth = "none", ...) {
  # --- 1. Input Validation ---
  if (!length(equate_recipe@design)) {
    cli::cli_abort("A design must be added with `add_design()` before adding a method.")
  }

  args <- list(...)

  method_norm <- tolower(method)
  smooth_norm <- tolower(smooth)
  design <- equate_recipe@design # "S", "R", or "CG"

  # --- 2. Normalize and Store Method ---
  new_method <- list()
  new_method$options <- args

  # Normalize primary method [cite: 135-138]
  new_method$method <- switch(method_norm,
                              "linear" = "L",
                              "equipercentile" = "E",
                              "irt" = "IRT",
                              cli::cli_abort("Unknown method: '{method}'. Choose 'linear', 'equipercentile', or 'irt'.")
  )

  # Normalize smoothing type [cite: 140-145]
  new_method$smooth <- switch(smooth_norm,
                              "none" = "N",
                              "beta_binomial" = "B",
                              "log_linear" = "L",
                              "cubic_spline" = "S",
                              "kernel" = "K",
                              "continuized_log_linear" = "Z",
                              cli::cli_abort("Unknown smooth type: '{smooth}'.")
  )

  if (new_method$smooth == "Z") {
    # This logic will be handled by get_method_options, but we can keep the alert
    if (is.null(args$cu) && is.null(args$cv) && is.null(args$cuv)) {
      cli::cli_alert_info("Using default log-linear smoothing parameters (continuized_log_linear). See documentation for details.")
    }
  }

  # Normalize calculation type based on method and design
  if (new_method$method == "L") {
    # The mean_only flag is now stored in options and handled by get_method_options
    if (design %in% c("R", "S")) {
      # The type is determined by mean_only, but we can store a user-friendly name
      is_mean <- args$mean_only %||% FALSE
      new_method$type <- ifelse(is_mean, "mean", "linear")
    } else { # design == "CG"
      new_method$type <- if (type == "default") "all" else tolower(type)
    }
  } else if (new_method$method == "E") {
    if (design %in% c("R", "S")) {
      new_method$type <- "identity" # No subtypes for RG/SG
    } else { # design == "CG"
      new_method$type <- switch(tolower(type),
                                "default" = "E", "e" = "E", "frequency" = "E",
                                "chained" = "C", "c" = "C",
                                "modified_frequency" = "F", "f" = "F",
                                "all" = "A", "g" = "G", "h" = "H",
                                cli::cli_abort("Unknown equipercentile type for CG design: '{type}'.")
      )
    }
  } else {
    new_method$type <- tolower(type)
  }

  # --- 3. Combination Validation ---
  if (new_method$smooth == "B" && design != "R") {
    cli::cli_abort("Beta-binomial smoothing is only supported for the Random Groups (RG) design.")
  }
  if (new_method$method == "L" && new_method$smooth != "N") {
    cli::cli_abort("Smoothing is not applicable to the linear method in Equating Recipes.")
  }

  # Create a unique title for the method
  type_str <- paste(new_method$type, collapse = "_")
  title_parts <- c(design, new_method$method, type_str, new_method$smooth)
  if (isTRUE(args$mean_only)) {
    title_parts <- c(title_parts, "mean_only")
  }
  new_method$title <- paste(title_parts, collapse = " ")

  if(new_method$title %in% names(equate_recipe@methods)) {
    cli::cli_alert_warning("A method with the same configuration already exists and will be overwritten.")
  }
  equate_recipe@methods[[new_method$title]] <- new_method

  equate_recipe
}


run_equating <- function(eq, boot_type = "perc", boot_replications = 1000){

  results <- lapply(1:nrow(eq@plan) |> `names<-`(apply(eq@plan, 1, paste0,collapse = ";")), \(i) {

    lapply(eq@methods, \(method){
      # Correctly extract the forms for the current iteration
      forms <- as.character(eq@plan[i,])

      # Extract method details from the unique title
      method_details <- strsplit(method$title, split = " ")[[1]]
      design <- method_details[1]
      method_code <- method_details[2]

      equate(forms = forms,
             design = design,
             type = method$type, # Use the stored type directly
             eq = eq,
             title = method$title,
             boot_type = boot_type,
             boot_replications = boot_replications,
             method = method_code)
    })
  })

  # Assign the custom class for S3 dispatch
  class(results) <- "equate_results"
  return(results)
}
