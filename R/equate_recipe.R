# Suggested S7 object definition
equate_recipe <- S7::new_class(
  "equate_recipe",
  properties = list(
    forms    = S7::class_list,
    data     = S7::class_list,
    data_ids = S7::class_list,
    plan     = S7::class_data.frame,
    design   = S7::class_character,
    methods  = S7::class_list,
    results = S7::class_list # We are actually going to store results here. We will then go to plot, print, and summary methods. Is this too crazy? I'm not sure how this will work out.
  )
)

.onLoad <- function(...) {
  S7::methods_register()
}

#' Constructor to initialize the recipe
#' @export
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
#' @param rel Optional. The internal consistency reliability (e.g., Cronbach's Alpha)
#'   of the form. If `NULL`, it is calculated automatically.
#'
#' @return The updated `equate_recipe` object containing the new form.
#' @export
add_form <- function(equate_recipe, crosstabs, name, id_cols = NULL, min_score = NULL, max_score = NULL, inc = NULL, rel = NULL) {

  data_cols <- !colnames(crosstabs) %in% id_cols
  form_data <- crosstabs[, data_cols]

  # Get the score scale attributes
  attrs <- get_form_attributes(form_data, form_name = name, min_score = min_score, max_score = max_score, inc = inc, rel = rel)

  # Safely attach the attributes to the data frame
  attr(form_data, "min") <- attrs$min
  attr(form_data, "max") <- attrs$max
  attr(form_data, "inc") <- attrs$inc
  attr(form_data, "range") <- attrs$range
  attr(form_data, "points") <- attrs$points
  attr(form_data, "n") <- attrs$n
  attr(form_data, "k") <- attrs$k
  attr(form_data, "rel") <- attrs$rel

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
#' \dontrun{
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
#' }
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
#' @param counter_balance_by Character vector of length 1. The name of the id_col to use in single group equating.
#' @return The updated `equate_recipe` object with the normalized `@design` field.
#' @export
#'
#' @examples
#' \dontrun{
#' equate_recipe |> add_design("single")
#' equate_recipe |> add_design("random-groups")
#' equate_recipe |> add_design("common-item nonequivalent")
#' }
add_design <- function(equate_recipe, design, counter_balance_by = NULL) {
  design <- tolower(design)

  # Match normalized design codes
  if (design %in% c("single", "single-group","single group", "s", "sg")) {
    cli::cli_inform("Single Group Design {ifelse(is.null(counter_balance_by), 'Without', 'With')} Counter Balancing")
    equate_recipe@design <- ("S") |> `attr<-`("counter_balance_by", counter_balance_by) |> `attr<-`("label", glue::glue("Single Group Design {ifelse(is.null(counter_balance_by), 'Without', 'With')} Counter Balancing"))
  } else if (design %in% c("random", "random-groups", "r", "rg")) {
    cli::cli_inform("Random/Equivalent Groups Design")
    equate_recipe@design <- "R"  |> `attr<-`("label", "Random/Equivalent Groups Design")
  } else if (design %in% c(
    "common-item", "common", "common-item nonequivalent", "common-item non-equivalent",
    "non-equivalent", "nonequivalent", "ne", "cg", "c", "n", "cneg"
  )) {
    cli::cli_inform("Common Item/Non-Equivalent Groups Design")
    equate_recipe@design <- "CG"  |> `attr<-`("label", "Common Item/Non-Equivalent Groups Design")
  } else {
    cli::cli_abort("Unknown design type: '{design}'")
  }

  equate_recipe
}

#' Get Master List of Default Method Options
#'
#' @description
#' Provides a master list of default arguments for all equating methods and their
#' options. This function serves as a single source of truth for default values
#' used throughout the package, primarily by `add_method()`.
#'
#' @return A list containing the complete set of default options.
#' @keywords internal
get_method_options <- function() {
  list(
    # General Bootstrap
    bootstrap = FALSE,
    boot_reps = 1000,
    boot_type = "perc",

    # General Linear
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
    rel = NULL,

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
    inc_v = 1,

    # IRT Options for all designs
    irt_pars = NULL, # Can be a single object or a named list for CNEG
    common_items = NULL, # For CNEG design
    transform_method = "stocking_lord", # For CNEG design
    anchor_type = "true_score", # For CNEG diagnostic
    theta = seq(-4, 4, length.out = 100),
    theta_dist = "normal", # Can be a single value or a named list for CNEG
    w1_irt = 0.5 # Default for IRT CNEG design
  )
}

#' Add an equating method to a recipe
#'
#' This function adds a specific equating method, including its calculation type
#' and smoothing options, to an `equate_recipe` object. You can add multiple
#' methods to a single recipe, and they will all be executed by `run_equating()`.
#'
#' @param equate_recipe An object of class `equate_recipe`.
#' @param method A character string specifying the primary equating method.
#'   Can be `"linear"`, `"equipercentile"`, `"irt"`, or `"identity"`.
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
#' \strong{`method = "identity"`}
#' \itemize{
#'   \item A baseline "equating" where the equated score is simply the original score.
#'   The `type` and `smooth` arguments are ignored.
#' }
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
#' \strong{`method = "irt"`}
#' \itemize{
#'   \item `type = "true_score"` (Default): Performs IRT true score equating.
#'   \item `type = "observed_score"`: Performs IRT observed score equating.
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
#'   \itemize{
#'     \item \code{nparm}: The number of parameters for the beta distribution (2 or 4). Integer. *Default: \code{4}*.
#'   }
#'   \item \code{smooth = "log_linear"}: Univariate log-linear presmoothing.
#'   \itemize{
#'     \item \code{degree}: The degree of the polynomial to fit. Integer. *Default: \code{3}*.
#'   }
#'   \item \code{smooth = "cubic_spline"}: Post-smoothing using cubic splines.
#'   \itemize{
#'     \item \code{s}: The smoothing parameter. Numeric. *Default: \code{0.2}*.
#'   }
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
#' \strong{IRT Method Options}
#' \itemize{
#'   \item \code{irt_pars}: This is a **required** argument.
#'     \itemize{
#'       \item For SG, RG, or concurrently calibrated CG designs, this should be a single `mirt` model object or a data frame containing parameters for all items.
#'       \item For the CNEG design with separately calibrated forms, this must be a **named list** of `mirt` objects or data frames, where the names match the form names in the recipe (e.g., `list(formX = pars_x, formY = pars_y)`).
#'     }
#'   \item \code{common_items}: For the CNEG design, a two-column data frame mapping the item identifiers for the common items. The column names must match the form names (e.g., `data.frame(formX = c("i1", "i3"), formY = c("i2", "i5"))`). If not provided, the function will default to using items with identical names on both forms.
#'   \item \code{transform_method}: For CNEG design, the method used to link the two ability scales. Can be `"stocking_lord"`, `"haebara"`, `"mean_sigma"`, or `"mean_mean"`. *Default: \code{"stocking_lord"}*.
#'   \item \code{anchor_type}: For the CNEG design, the equating type to use for the diagnostic equating of the anchor test. If not specified, this will default to the main `type` (e.g., `"true_score"` or `"observed_score"`).
#'   \item \code{theta}: A numeric vector defining the grid for the latent ability scale. *Default: \code{seq(-4, 4, length.out = 100)}*.
#'   \item \code{theta_dist}: For observed score equating.
#'     \itemize{
#'       \item For SG/RG designs, a character string (`"normal"` or `"uniform"`) or a numeric vector defining the single population ability distribution. *Default: \code{"normal"}*.
#'       \item For CNEG design, this must be a **named list** where names match the form names and values are the distributions for each group (e.g., `list(formX = "normal", formY = "normal")`).
#'     }
#'   \item \code{w1}: For CNEG observed score equating, the weight given to the new group (Form X) in the synthetic population. *Default: \code{0.5}*.
#'   \item \code{bootstrap}: Should parametric bootstrap standard errors be calculated? Logical. *Default: \code{FALSE}*. Requires a `mirt` model object for `irt_pars`.
#'   \item \code{boot_reps}: The number of bootstrap replications. Integer. *Default: \code{100}* for IRT, \code{1000} otherwise.
#' }
#'
#' @examples
#' \dontrun{
#' recipe <- init_equating() |> add_design("common-item")
#' # Add a linear method that runs both Tucker and Chained equating
#' recipe <- recipe |> add_method("linear", type = c("tucker", "chained"))
#' # Add an identity equating as a baseline
#' recipe <- recipe |> add_method("identity")
#' }
add_method <- function(equate_recipe, method, type = "default", smooth = "none", ...) {
  # --- 1. Setup ---
  if (!length(equate_recipe@design)) {
    cli::cli_abort("A design must be added with `add_design()` before adding a method.")
  }

  user_args <- list(...)
  method_norm <- tolower(method)
  smooth_norm <- tolower(smooth)
  design <- equate_recipe@design

  # --- 2. Consolidate and Validate Options ---

  # Get master list of defaults
  defaults <- get_method_options()

  # Apply special-case defaults based on method or smoothing
  if (method_norm == "irt") {
    defaults$boot_reps <- 100
    # Use the IRT-specific w1 default
    defaults$w1 <- defaults$w1_irt
  }
  if (smooth_norm == "continuized_log_linear" && is.null(user_args$boot_reps)) {
    defaults$boot_reps <- 1
  }

  # Merge user arguments with defaults so user args take precedence
  options <- defaults
  for (name in names(user_args)) {
    options[[name]] <- user_args[[name]]
  }

  new_method <- list()
  new_method$method <- switch(method_norm,
                              "linear" = "L",
                              "equipercentile" = "E",
                              "irt" = "IRT",
                              "identity" = "I",
                              cli::cli_abort("Unknown method: '{method}'. Choose 'linear', 'equipercentile', 'irt', or 'identity'.")
  )

  # Normalize calculation type based on method and design
  if (new_method$method == "I") {
    new_method$type <- "identity"
  } else if (new_method$method == "L") {
    # if(new_method$design == "SG") new_method$type <- ifelse(options$mean_only, "mean", "linear") #### Gemini, this is a poor design choice.
    new_method$type <- "all"


  } else if (new_method$method == "E") {
    new_method$type <- "E" # Placeholder
  } else if (new_method$method == "IRT") {
    new_method$type <- switch(tolower(type),
                              "default" = "true_score",
                              "true_score" = "true_score",
                              "observed_score" = "observed_score",
                              cli::cli_abort("Unknown IRT type: '{type}'. Choose 'true_score' or 'observed_score'.")
    )
    # If anchor_type was not specified by user, make it default to the main IRT type
    if (is.null(user_args$anchor_type)) {
      options$anchor_type <- new_method$type
    }
  } else {
    new_method$type <- tolower(type)
  }

  # --- Method-Specific Validations ---
  if (new_method$method == "IRT") {
    if (is.null(options$irt_pars)) {
      cli::cli_abort("For the 'irt' method, the 'irt_pars' argument must be provided.")
    }

    if (options$bootstrap) {
      is_mirt_object <- function(obj) (attr(class(obj), "package") %||% "not") == "mirt"
      pars_to_check <- options$irt_pars
      if (is.list(pars_to_check) && !is.data.frame(pars_to_check)) {
        if (!all(sapply(pars_to_check, is_mirt_object))) {
          cli::cli_abort("For bootstrapped IRT standard errors, parameters must be in 'mirt' model objects.")
        }
      } else if (!is_mirt_object(pars_to_check)) {
        cli::cli_abort("For bootstrapped IRT standard errors, parameters must be in a 'mirt' model object.")
      }
      cli::cli_alert_info("Parametric bootstrapping for IRT methods is computationally intensive and may take a long time.")
    }

    if (!is.numeric(options$theta)) {
      cli::cli_abort("'theta' must be a numeric vector.")
    }

    # Store only the relevant options for this method
    final_options <- options[c(
      "irt_pars", "theta", "theta_dist", "w1", "common_items",
      "transform_method", "anchor_type", "bootstrap", "boot_reps"
    )]
  } else {
    final_options <- options
  }

  # Store the final, complete set of options
  new_method$options <- final_options


  # Normalize smoothing type
  new_method$smooth <- switch(smooth_norm,
                              "none" = "N",
                              "beta_binomial" = "B",
                              "log_linear" = "L",
                              "cubic_spline" = "S",
                              "kernel" = "K",
                              "continuized_log_linear" = "Z",
                              cli::cli_abort("Unknown smooth type: '{smooth}'.")
  )


  # --- 3. Combination Validation ---
  if (new_method$method %in% c("IRT", "I") && new_method$smooth != "N") {
    cli::cli_abort("Smoothing is not applicable to the '{method}' method.")
  }
  if (new_method$method == "L" && new_method$smooth != "N") {
    cli::cli_abort("Smoothing is not applicable to the linear method.")
  }

  # Create a unique title for the method
  type_str <- paste(new_method$type, collapse = "_")
  title_parts <- c(design, new_method$method, type_str, new_method$smooth)
  new_method$title <- paste(title_parts, collapse = " ")

  if (new_method$title %in% names(equate_recipe@methods)) {
    cli::cli_alert_warning("A method with the same configuration already exists and will be overwritten.")
  }

  equate_recipe@methods[[new_method$title]] <- new_method

  equate_recipe
}







run_equating <- function(eq){#, boot_type = "perc", boot_replications = 1000

  eq@results <- lapply(1:nrow(eq@plan) |> `names<-`(apply(eq@plan, 1, paste0,collapse = ";")), \(i) {

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
             method = method_code)
    })
  })

  # Assign the custom class for S3 dispatch
  print(summary(eq)) # So we need to make a method for summary for this object.
  # It will need to
  return(eq)
}
