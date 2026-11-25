# Sets of specs: the workflowsets analogue ---------------------------------------------------------

#' A set of simulation specs to compare
#'
#' Like `workflowsets::workflow_set()`. Useful when candidates take different
#' parameters, for example RIL designs (tune `selfing`) versus NAM designs
#' (tune `families`) versus MAGIC (tune `founders`): each spec gets its own
#' grid, but all run on the same replicate seeds and are ranked together.
#'
#' @param ... Named simulation specs.
#' @return A tibble of class `sim_set` with columns `wflow_id`, `spec` and `result`.
#' @examples
#' sim_set(
#'   ril = mapping_population(type = "RIL", lines = tune()),
#'   nam = mapping_population(type = "NAM", lines = tune(), families = tune())
#' )
#' @export
sim_set <- function(...) {
  specs <- list(...)
  if (!length(specs) || is.null(names(specs)) || any(names(specs) == "")) abort("Name every spec.")
  if (!all(vapply(specs, inherits, logical(1), "sim_spec"))) abort("Every element must be a simulation spec.")
  out <- tibble::tibble(wflow_id = names(specs), spec = unname(specs), result = rep(list(NULL), length(specs)))
  structure(out, class = c("sim_set", class(out)))
}

#' Run a tuner over every spec in a set
#'
#' Like `workflowsets::workflow_map()`. Every spec uses the same replicates,
#' so designs from different specs are paired by seed.
#'
#' @param object A [sim_set()].
#' @param fn Tuner name: `"tune_grid"`, `"tune_race"`, `"tune_bayes"` or `"fit_resamples"`.
#' @param ... Arguments passed to the tuner for every spec (e.g. `resamples`, `metrics`, `grid`).
#' @param grids Optional named list of grids, one per `wflow_id`, overriding `grid`.
#' @param verbose Print which spec is running.
#' @return The set with its `result` column filled.
#' @export
sim_map <- function(object, fn = "tune_grid", ..., grids = NULL, verbose = FALSE) {
  if (!inherits(object, "sim_set")) abort("`object` must come from sim_set().")
  fn <- rlang::arg_match0(fn, c("tune_grid", "tune_race", "tune_bayes", "fit_resamples"))
  tuner <- switch(fn, tune_grid = tune_grid.sim_spec, tune_race = tune_race, tune_bayes = tune_bayes.sim_spec,
                  fit_resamples = fit_resamples.sim_spec)
  dots <- list(...)
  for (i in seq_len(nrow(object))) {
    id <- object$wflow_id[i]
    if (verbose) cli::cli_inform("Running {fn} for {id}.")
    args <- dots
    if (!is.null(grids[[id]])) args$grid <- grids[[id]]
    spec <- object$spec[[i]]
    if (fn != "fit_resamples" && !nrow(tune_args.sim_spec(spec))) {
      args$grid <- NULL
      args$param_info <- NULL
      object$result[[i]] <- do.call(fit_resamples.sim_spec,
                                    c(list(spec), args[intersect(names(args), c("resamples", "metrics", "control"))]))
    } else {
      object$result[[i]] <- do.call(tuner, c(list(spec), args))
    }
  }
  object
}

#' @rdname collect_metrics
#' @export
collect_metrics.sim_set <- function(x, ..., summarize = TRUE, by_scenario = NULL, type = c("long", "wide")) {
  type <- match.arg(type)
  pieces <- lapply(seq_len(nrow(x)), function(i) {
    res <- x$result[[i]]
    if (is.null(res)) abort(sprintf("`%s` has not been run; call sim_map() first.", x$wflow_id[i]))
    m <- collect_metrics.sim_results(res, summarize = summarize, by_scenario = by_scenario, type = type)
    ids <- res$params$id
    params <- if (length(ids)) {
      vapply(seq_len(nrow(m)), function(r) paste(ids, vapply(ids, function(id) format(m[[id]][r]), ""), sep = "=",
                                                 collapse = ", "), "")
    } else rep("", nrow(m))
    keep <- setdiff(names(m), ids)
    tibble::as_tibble(cbind(tibble::tibble(wflow_id = x$wflow_id[i], parameters = params), m[, keep, drop = FALSE]))
  })
  do.call(rbind, pieces)
}

#' Rank designs across a set of specs
#'
#' Like `workflowsets::rank_results()`.
#'
#' @param x A [sim_set()] that has been run with [sim_map()].
#' @param rank_metric Metric to rank by; defaults to the first in the set.
#' @param constraints Named list of bounds on metric means, as in [show_best()].
#' @param select_best Keep only the best design of each spec.
#' @return A wide tibble, best first, with a `rank` column.
#' @export
rank_results <- function(x, rank_metric = NULL, constraints = NULL, select_best = FALSE) {
  if (!inherits(x, "sim_set")) abort("`x` must come from sim_set().")
  wide <- collect_metrics.sim_set(x, by_scenario = FALSE, type = "wide")
  first <- x$result[[1]]
  metrics <- first$metrics
  m <- if (is.null(rank_metric)) metrics[[1]] else metrics[[rank_metric]]
  if (is.null(m)) abort(sprintf("`%s` is not in the metric set.", rank_metric))
  ok <- rep(TRUE, nrow(wide))
  for (name in names(constraints)) {
    rule <- constraint_rule(constraints[[name]], metrics[[name]])
    v <- wide[[name]]
    ok <- ok & is.finite(v) & switch(rule$op, ">=" = v >= rule$bound, ">" = v > rule$bound,
                                     "<=" = v <= rule$bound, "<" = v < rule$bound)
  }
  wide <- wide[ok, , drop = FALSE]
  key <- if (identical(m$direction, "zero")) abs(wide[[m$name]]) else -metric_sign(m) * wide[[m$name]]
  wide <- wide[order(key), , drop = FALSE]
  if (select_best) wide <- wide[!duplicated(wide$wflow_id), , drop = FALSE]
  wide$rank <- seq_len(nrow(wide))
  wide[, c("wflow_id", "rank", setdiff(names(wide), c("wflow_id", "rank")))]
}

#' Pull one spec's results out of a set
#' @param x A run [sim_set()].
#' @param id A `wflow_id`.
#' @return A `sim_results` object.
#' @export
extract_sim_result <- function(x, id) {
  i <- match(id, x$wflow_id)
  if (is.na(i)) abort(sprintf("No spec `%s` in the set.", id))
  x$result[[i]]
}

#' @export
print.sim_set <- function(x, ...) {
  cat(sprintf("# A simulation set of %d spec(s)\n", nrow(x)))
  run <- !vapply(x$result, is.null, logical(1))
  info <- tibble::tibble(wflow_id = x$wflow_id, spec = vapply(x$spec, function(s) class(s)[1], ""),
                         tuned = vapply(x$spec, function(s) paste(tune_args.sim_spec(s)$id, collapse = ", "), ""),
                         result = ifelse(run, vapply(x$result, function(r) if (is.null(r)) "" else
                           sprintf("%s: %d designs", r$method, nrow(r$configs)), ""), "<not run>"))
  print(info)
  invisible(x)
}
