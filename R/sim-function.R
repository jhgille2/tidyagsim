#' A simulation spec from any R function
#'
#' Wraps a function `fn(..., seed)` that runs one replicate and returns a named
#' list (or one-row data frame) of outcomes, or a single number (stored as
#' `value`). Arguments given here are passed to `fn`; set any of them to
#' `tune()` and give its range in `params`.
#'
#' @param fn A function whose arguments include the ones given in `...` and `seed`.
#' @param ... Arguments for `fn`; use `tune()` to mark one for tuning.
#' @param params Named list of dials parameter objects giving default ranges for
#'   tuned arguments (see [sim_param()]).
#' @param check Optional function of the same arguments returning `NULL` when a
#'   design is feasible or a reason string when it is not.
#' @param metrics Optional default [sim_metric_set()].
#' @return A `sim_function` spec.
#' @examples
#' toy <- sim_function(
#'   function(x, noise, seed) {
#'     set.seed(seed)
#'     list(y = -(x - 3)^2 + rnorm(1, sd = noise))
#'   },
#'   x = tune(), noise = 0.5,
#'   params = list(x = sim_param(c(0, 6)))
#' )
#' res <- tune_grid(toy, sim_replicates(5), grid = 7,
#'                  metrics = sim_metric_set(sim_metric("y", "maximize")))
#' show_best(res, n = 2)
#' @export
sim_function <- function(fn, ..., params = list(), check = NULL, metrics = NULL) {
  if (!is.function(fn)) abort("`fn` must be a function.")
  if (!"seed" %in% names(formals(fn)) && !"..." %in% names(formals(fn))) {
    abort("`fn` must take a `seed` argument.")
  }
  new_sim_spec("sim_function", list(...), engine = "R", fn = fn, params = params, check = check,
               default_metrics = metrics)
}

#' @export
sim_arg_names.sim_function <- function(object) {
  f <- names(formals(object$fn))
  if ("..." %in% f) NULL else setdiff(f, "seed")
}

#' @export
sim_parameters.sim_function <- function(object, ...) object$params

#' @export
sim_check.sim_function <- function(object, args, ...) {
  if (is.null(object$check)) return(NULL)
  do.call(object$check, args)
}

#' @export
sim_default_metrics.sim_function <- function(object) object$default_metrics

#' @export
sim_simulate.sim_function <- function(object, args, seeds, ...) {
  lapply(seeds, function(seed) {
    out <- with_seed(seed, do.call(object$fn, c(args, list(seed = seed))))
    if (is.data.frame(out)) out <- as.list(out[1, , drop = FALSE])
    if (is.numeric(out) && length(out) == 1 && is.null(names(out))) out <- list(value = out)
    as.list(out)
  })
}
