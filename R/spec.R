# Simulation specifications ------------------------------------------------------------

#' Simulation specifications
#'
#' A simulation spec is the tidyagsim analogue of a parsnip model specification.
#' It records the arguments of one simulator. Any argument may be set to
#' [tune()] to let a tuner choose it; the rest stay fixed. Every spec can be
#' passed to [tune_grid()], [tune_race()], [tune_bayes()] and [fit_resamples()].
#'
#' Spec constructors shipped with the package:
#'
#' * [mapping_population()]: QTL mapping populations, simulated in R.
#' * [sim_function()]: any R function of the arguments and a seed.
#' * [breeding_program()]: a breedsim flowchart program run by AlphaSimR.
#'
#' New simulators are added by writing a constructor that calls
#' `new_sim_spec()` and methods for [sim_simulate()] and, optionally,
#' [sim_check()] and [sim_parameters()].
#'
#' @param cls Subclass name.
#' @param args Named list of arguments; values may be `tune()` calls.
#' @param engine Engine name (the simulator implementation).
#' @param ... Other fields stored on the spec.
#' @return An object of class `c(cls, "sim_spec")`.
#' @name sim_spec
#' @export
new_sim_spec <- function(cls, args, engine = "tidyagsim", ...) {
  if (!is.list(args) || (length(args) && (is.null(names(args)) || any(names(args) == "")))) {
    abort("`args` must be a named list.")
  }
  structure(list(args = args, engine = engine, ...), class = c(cls, "sim_spec"))
}

#' Is an argument a tune() placeholder?
#' @param x An object.
#' @return `TRUE` or `FALSE`.
#' @export
is_tune <- function(x) {
  is.call(x) && identical(x[[1]], quote(tune))
}

tune_id <- function(x, name) {
  if (length(x) > 1 && is.character(x[[2]]) && nzchar(x[[2]])) x[[2]] else name
}

#' Tuned arguments of a spec
#'
#' @param object A simulation spec.
#' @param ... Not used.
#' @return A tibble with one row per argument marked `tune()`: `name`, `tunable`, `id`,
#'   `source` and `component`, as in tidymodels' `tune_args()`.
#' @exportS3Method generics::tune_args
tune_args.sim_spec <- function(object, ...) {
  marked <- Filter(is_tune, object$args)
  tibble::tibble(
    name = names(marked),
    tunable = rep(TRUE, length(marked)),
    id = vapply(names(marked), function(n) tune_id(marked[[n]], n), character(1), USE.NAMES = FALSE),
    source = rep("sim_spec", length(marked)),
    component = rep(class(object)[1], length(marked)),
    component_id = rep("main", length(marked))
  )
}

#' Default parameter objects of a spec
#'
#' Returns a named list of dials parameter objects, one per argument that has a
#' sensible default range. Spec authors implement this generic so that
#' `tune()` works without the user writing ranges.
#'
#' @param object A simulation spec.
#' @param ... Not used.
#' @return A named list of `param` objects.
#' @export
sim_parameters <- function(object, ...) UseMethod("sim_parameters")

#' @export
sim_parameters.default <- function(object, ...) list()

#' Parameters that can be tuned, like tidymodels' tunable()
#'
#' @param x A simulation spec.
#' @param ... Not used.
#' @return A tibble with `name`, `call_info`, `source`, `component` and `component_id`.
#' @exportS3Method generics::tunable
tunable.sim_spec <- function(x, ...) {
  defaults <- sim_parameters(x)
  tibble::tibble(
    name = names(defaults),
    call_info = lapply(names(defaults), function(n) list(pkg = "tidyagsim", fun = n)),
    source = "sim_spec",
    component = class(x)[1],
    component_id = "main"
  )
}

#' Extract the parameter set of a spec
#'
#' The parameter set of the arguments marked [tune()], with default ranges from
#' the spec. Change ranges with `update()`, exactly as in tidymodels, and pass
#' the result to a tuner as `param_info`.
#'
#' @param x A simulation spec.
#' @param ... Not used.
#' @return A dials `parameters` object.
#' @examples
#' spec <- mapping_population(type = "RIL", lines = tune(), selfing = tune())
#' params <- extract_parameter_set_dials(spec)
#' params
#' update(params, lines = n_lines(c(100, 400)))
#' @exportS3Method hardhat::extract_parameter_set_dials
extract_parameter_set_dials.sim_spec <- function(x, ...) {
  marked <- tune_args.sim_spec(x)
  marked <- marked[!duplicated(marked$id), , drop = FALSE]
  if (!nrow(marked)) {
    return(dials::parameters(list()))
  }
  defaults <- sim_parameters(x)
  objects <- lapply(seq_len(nrow(marked)), function(i) {
    name <- marked$name[i]
    obj <- defaults[[name]]
    if (is.null(obj)) {
      abort(c(sprintf("No default range for `%s`.", name),
              i = "Supply one with `update(extract_parameter_set_dials(spec), ...)` and pass it as `param_info`."))
    }
    obj
  })
  dials::parameters_constr(
    name = marked$name, id = marked$id, source = marked$source,
    component = marked$component, component_id = marked$component_id,
    object = objects
  )
}

#' Change the arguments of a spec
#'
#' @param object A simulation spec.
#' @param ... Named arguments; use `tune()` to mark one for tuning.
#' @return The modified spec.
#' @examples
#' mapping_population(type = "F2") |> set_args(lines = tune())
#' @export
set_args <- function(object, ...) UseMethod("set_args")

#' @export
set_args.default <- function(object, ...) {
  forward("parsnip", "set_args", object, ...)
}

#' @export
set_args.sim_spec <- function(object, ...) {
  new <- list(...)
  if (length(new) && (is.null(names(new)) || any(names(new) == ""))) abort("Name every argument.")
  check_arg_names(object, names(new))
  object$args[names(new)] <- new
  object
}

#' Change the engine of a spec
#' @param object A simulation spec.
#' @param engine Engine name.
#' @param ... Engine-specific options stored on the spec.
#' @return The modified spec.
#' @export
set_engine <- function(object, engine, ...) UseMethod("set_engine")

#' @export
set_engine.default <- function(object, engine, ...) forward("parsnip", "set_engine", object, engine, ...)

#' @export
set_engine.sim_spec <- function(object, engine, ...) {
  object$engine <- engine
  object$engine_args <- utils::modifyList(object$engine_args %||% list(), list(...))
  object
}

#' Argument names a spec accepts
#' @param object A simulation spec.
#' @return A character vector, or `NULL` when any name is accepted.
#' @keywords internal
#' @export
sim_arg_names <- function(object) UseMethod("sim_arg_names")

#' @export
sim_arg_names.default <- function(object) NULL

check_arg_names <- function(object, names) {
  known <- sim_arg_names(object)
  bad <- setdiff(names, known)
  if (!is.null(known) && length(bad)) {
    abort(sprintf("%s has no argument%s %s. Known: %s.", class(object)[1], if (length(bad) > 1) "s" else "",
                  paste0("`", bad, "`", collapse = ", "), paste(known, collapse = ", ")))
  }
  invisible(TRUE)
}

#' Concrete arguments for one design
#'
#' Replaces each `tune()` argument by the design's value (looked up by id) and
#' applies scenario overrides.
#'
#' @param object A simulation spec.
#' @param values Named list of tuned values, keyed by parameter id.
#' @param scenario Named list of scenario values (may be empty).
#' @return A named list of arguments.
#' @keywords internal
#' @export
resolve_args <- function(object, values = list(), scenario = list()) {
  marked <- tune_args.sim_spec(object)
  args <- object$args
  missing <- setdiff(marked$id, names(values))
  if (length(missing)) abort(sprintf("No value for tuned parameter(s) %s.", paste(missing, collapse = ", ")))
  for (i in seq_len(nrow(marked))) args[[marked$name[i]]] <- values[[marked$id[i]]]
  if (length(scenario)) {
    check_arg_names(object, names(scenario))
    args[names(scenario)] <- scenario
  }
  args
}

#' Finalize a spec with chosen values
#'
#' Like `finalize_model()`: replaces every `tune()` argument with the value in
#' `parameters` (a one-row tibble, for example from [select_best()]).
#'
#' @param x A simulation spec.
#' @param parameters A one-row data frame or named list of values keyed by id.
#' @return A spec with no tuned arguments.
#' @export
finalize_sim <- function(x, parameters) {
  if (!inherits(x, "sim_spec")) abort("`x` must be a simulation spec.")
  values <- as.list(parameters)
  values <- lapply(values, function(v) if (is.factor(v)) as.character(v) else v)
  x$args <- resolve_args(x, values)
  x
}

#' Run one design of a spec
#'
#' The generic every simulator implements. It receives fully resolved arguments
#' and a vector of seeds and returns a list with one named list of outcomes per
#' seed (for example `list(power = 1, cost = 9800)`).
#'
#' @param object A simulation spec.
#' @param args Resolved arguments (no `tune()` left).
#' @param seeds Integer seeds, one per replicate.
#' @param ... Not used.
#' @return A list of named lists, one per seed.
#' @export
sim_simulate <- function(object, args, seeds, ...) UseMethod("sim_simulate")

#' Check whether a design is feasible
#'
#' @param object A simulation spec.
#' @param args Resolved arguments.
#' @param ... Not used.
#' @return `NULL` if the design can be run, otherwise a one-line reason.
#' @export
sim_check <- function(object, args, ...) UseMethod("sim_check")

#' @export
sim_check.default <- function(object, args, ...) NULL

#' Default metric set of a spec
#' @param object A simulation spec.
#' @return A metric set, or `NULL`.
#' @keywords internal
#' @export
sim_default_metrics <- function(object) UseMethod("sim_default_metrics")

#' @export
sim_default_metrics.default <- function(object) NULL

#' @export
print.sim_spec <- function(x, ...) {
  cat(sprintf("%s simulation spec (engine: %s)\n", class(x)[1], x$engine))
  marked <- vapply(x$args, is_tune, logical(1))
  if (any(marked)) {
    cat("\nTuned arguments:\n")
    for (n in names(x$args)[marked]) cat(sprintf("  %s = %s\n", n, deparse(x$args[[n]])))
  }
  if (any(!marked)) {
    cat("\nFixed arguments:\n")
    for (n in names(x$args)[!marked]) cat(sprintf("  %s = %s\n", n, format_arg(x$args[[n]])))
  }
  invisible(x)
}

format_arg <- function(v) {
  if (is.function(v)) return("<function>")
  if (is.list(v)) return(sprintf("<list of %d>", length(v)))
  if (length(v) > 6) return(paste0(paste(format(utils::head(v, 6)), collapse = ", "), ", ..."))
  paste(format(v), collapse = ", ")
}

forward <- function(pkg, fun, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort(sprintf("No tidyagsim method for this object, and package %s is not installed.", pkg))
  }
  getExportedValue(pkg, fun)(...)
}

spec_fingerprint <- function(object) {
  clean <- object
  clean$args <- lapply(object$args, function(a) if (is.function(a)) deparse(a) else a)
  if (!is.null(clean$fn)) clean$fn <- deparse(clean$fn)
  rlang::hash(clean)
}
