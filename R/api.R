# ---------------------------------------------------------------------------
# ModelsCloud API surface for BCRAT.
#
#   - model_run(model_input)  : synchronous; 5-year risk for each input row
#   - get_sample_input(n)     : an example cohort
#   - get_default_input()     : one baseline woman to modify
# ---------------------------------------------------------------------------

# Predictor columns accepted by the API. hispanic_nativity and asian_subgroup
# are optional (only meaningful for Hispanic/Asian `race`, and default
# sensibly in bcrat() when absent); atypical_hyperplasia also defaults.
.bcrat_required_vars <- c(
  "age", "race", "age_menarche", "age_first_birth", "n_biopsies", "n_relatives"
)
.bcrat_optional_vars <- c("atypical_hyperplasia", "hispanic_nativity", "asian_subgroup")
.bcrat_vars <- c(.bcrat_required_vars, .bcrat_optional_vars)

# Accepted aliases (alias -> canonical), including NCI's own raw field names.
.bcrat_alias <- c(
  t1 = "age",
  agemen = "age_menarche", menarche = "age_menarche",
  age1st = "age_first_birth", first_birth = "age_first_birth",
  age_at_first_birth = "age_first_birth",
  n_biop = "n_biopsies", num_biopsies = "n_biopsies", biopsies = "n_biopsies",
  n_rels = "n_relatives", num_relatives = "n_relatives", relatives = "n_relatives",
  family_hist_count = "n_relatives",
  hypplas = "atypical_hyperplasia", hyperplasia = "atypical_hyperplasia",
  nativity = "hispanic_nativity",
  ethnicity = "asian_subgroup", subgroup = "asian_subgroup"
)

# Normalise an incoming payload: accept either a wrapped `model_input` object
# or fields passed directly as named arguments (`dots`); rename known aliases
# to canonical names; drop unrecognised extra fields. Returns NULL when no
# input at all was supplied.
.bcrat_normalize <- function(model_input, dots) {
  if (is.null(model_input)) {
    if (length(dots) == 0) return(NULL)
    model_input <- dots
  }
  # A JSON `null` arrives as a length-0 NULL, which as.data.frame() cannot
  # reconcile with the length-1 fields beside it ("differing number of rows").
  # This matters because the model's own output round-trips: hispanic_nativity
  # and asian_subgroup are NA for anyone who isn't Hispanic/Asian, serialise as
  # null, and come straight back in on the next call. An explicit null and an
  # absent field mean the same thing here, so drop them and let the defaults in
  # bcrat() apply.
  if (is.list(model_input) && !is.data.frame(model_input)) {
    model_input <- model_input[!vapply(model_input, is.null, logical(1))]
    if (length(model_input) == 0) return(NULL)
  }
  df <- as.data.frame(model_input, stringsAsFactors = FALSE)
  names(df) <- tolower(names(df))
  for (a in intersect(names(df), names(.bcrat_alias))) {
    canon <- .bcrat_alias[[a]]
    if (!canon %in% names(df)) names(df)[match(a, names(df))] <- canon
  }
  df[, intersect(names(df), .bcrat_vars), drop = FALSE]
}

#' Run the BCRAT model (ModelsCloud entry point)
#'
#' Scores one or more women and returns the input augmented with the
#' predicted 5-year invasive breast cancer risk. This is the synchronous
#' prediction pattern expected by ModelsCloud: a table of women in, the same
#' table plus predictions out.
#'
#' @details
#' The function is deliberately forgiving about how inputs arrive, so it
#' works whether the platform wraps the fields under `model_input` or passes
#' them directly:
#'
#' * **Wrapped** — `model_run(model_input = list(age = 50, ...))` or a data
#'   frame (the form produced by the `modelscloud` client).
#' * **Unwrapped** — `model_run(age = 50, race = "white", ...)` (the form
#'   produced by a raw `do.call(model_run, funcInput)`).
#'
#' Common aliases (including NCI's own raw field names) are accepted and
#' mapped to the canonical names: `T1` -> `age`, `AgeMen` -> `age_menarche`,
#' `Age1st` -> `age_first_birth`, `N_Biop` -> `n_biopsies`, `N_Rels` ->
#' `n_relatives`, `HypPlas` -> `atypical_hyperplasia`, `nativity` ->
#' `hispanic_nativity`, `ethnicity`/`subgroup` -> `asian_subgroup`.
#' Unrecognised extra fields are ignored.
#'
#' @param model_input A named list (one woman) or data frame (one row per
#'   woman) whose columns are the BCRAT predictors. See [bcrat()] for the
#'   meaning and coding of each field, or call [get_sample_input()] /
#'   [get_default_input()] for ready-made examples. If `NULL` and no fields
#'   are supplied via `...`, the model's [get_default_input()] is used.
#' @param ... Alternative to `model_input`: the predictor fields supplied
#'   directly as named arguments (e.g. from an unwrapped API call).
#'
#' @return A data frame: the input columns plus `risk` (5-year probability in
#'   `[0, 1]`) and `risk_percent` (the same value as a percentage, rounded to
#'   two decimals).
#'
#' @seealso [bcrat()], [get_sample_input()], [get_default_input()]
#' @examples
#' model_run(get_sample_input())
#' model_run(get_default_input())
#' # Unwrapped + NCI-style aliases:
#' model_run(T1 = 50, race = "white", AgeMen = 13, Age1st = 24, N_Biop = 0,
#'           N_Rels = 0)
#' @export
model_run <- function(model_input = NULL, ...) {
  df <- .bcrat_normalize(model_input, list(...))
  if (is.null(df)) df <- as.data.frame(get_default_input(), stringsAsFactors = FALSE)

  missing <- setdiff(.bcrat_required_vars, names(df))
  if (length(missing) > 0) {
    stop("Missing required variable(s): ", paste(missing, collapse = ", "),
         ". Accepted names (incl. aliases) are documented in ?model_run.",
         call. = FALSE)
  }

  df$risk <- bcrat(
    age                   = df$age,
    race                  = df$race,
    age_menarche          = df$age_menarche,
    age_first_birth       = df$age_first_birth,
    n_biopsies            = df$n_biopsies,
    n_relatives           = df$n_relatives,
    atypical_hyperplasia  = if (is.null(df$atypical_hyperplasia)) "unknown" else df$atypical_hyperplasia,
    hispanic_nativity     = df$hispanic_nativity,
    asian_subgroup        = df$asian_subgroup
  )
  df$risk_percent <- round(100 * df$risk, 2)
  df
}

#' Example BCRAT input cohort
#'
#' Returns a small data frame of example women that can be passed straight to
#' [model_run()], covering each supported race/ethnicity category.
#'
#' @param n Optional positive integer; if supplied, the first `n` rows are
#'   returned. Defaults to all rows.
#' @param ... Additional fields supplied by the platform; ignored.
#' @return A data frame of example women with the BCRAT predictor columns.
#' @seealso [model_run()], [get_default_input()]
#' @examples
#' get_sample_input()
#' get_sample_input(n = 2)
#' @export
get_sample_input <- function(n = NULL, ...) {
  df <- data.frame(
    age                   = c(50, 45, 62, 58, 70),
    race                  = c("White", "Black", "Hispanic", "Asian", "American Indian or Alaska Native"),
    age_menarche          = c(13, 12, 14, 11, 13),
    age_first_birth       = c(24, "nulliparous", 26, 22, 20),
    n_biopsies            = c(0, 1, 0, 2, 1),
    n_relatives           = c(0, 1, 1, 0, 2),
    atypical_hyperplasia  = c("unknown", "no", "unknown", "yes", "no"),
    hispanic_nativity     = c(NA, NA, "foreign_born", NA, NA),
    asian_subgroup        = c(NA, NA, NA, "chinese", NA),
    stringsAsFactors      = FALSE
  )
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || n < 1L) {
      stop("`n` must be a single positive integer.", call. = FALSE)
    }
    df <- utils::head(df, n)
  }
  df
}

#' Default BCRAT input
#'
#' Returns a single baseline woman as a named list, ready to modify and pass
#' to [model_run()].
#'
#' @param ... Additional fields supplied by the platform; ignored.
#' @return A named list of default predictor values.
#' @seealso [model_run()], [get_sample_input()]
#' @examples
#' woman <- get_default_input()
#' woman$age <- 55
#' model_run(woman)
#' @export
get_default_input <- function(...) {
  list(
    age                   = 50,
    race                  = "White",
    age_menarche          = 13,
    age_first_birth       = 24,
    n_biopsies            = 0,
    n_relatives           = 0,
    atypical_hyperplasia  = "unknown"
  )
}
