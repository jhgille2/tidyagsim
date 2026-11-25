# Metric sets: the yardstick analogue ----------------------------------------------------

builtin_metrics <- function() {
  rows <- list(
    # Mapping populations
    c("power", "maximize", "QTL detection power", "proportion"),
    c("false_positive", "minimize", "Replicates with a false QTL peak", "proportion"),
    c("location_error", "minimize", "Distance from peak to true QTL", "cM"),
    c("support_interval", "minimize", "1.5-LOD support interval", "cM"),
    c("peak_lod", "maximize", "LOD score at the QTL", "LOD"),
    c("realized_pve", "maximize", "Realized variance explained by the QTL", "proportion"),
    c("threshold", "minimize", "Genome-wide LOD threshold", "LOD"),
    # Breeding programs
    c("release_gain", "maximize", "Genetic gain of released lines", "trait units"),
    c("gain_per_year", "maximize", "Genetic gain per year", "trait units per year"),
    c("gain_per_cost", "maximize", "Genetic gain per unit cost", "trait units per currency unit"),
    c("expected_heterozygosity", "maximize", "Parent-pool expected heterozygosity", "proportion"),
    c("genetic_variance", "maximize", "Additive genetic variance", "trait units squared"),
    # Resources
    c("lines", "minimize", "Lines", "lines"),
    c("plots", "minimize", "Field plots", "plots"),
    c("genotyped", "minimize", "Lines genotyped", "lines"),
    c("crosses", "minimize", "Crosses made", "crosses"),
    c("cost", "minimize", "Total cost", "currency"),
    c("years", "minimize", "Years from first cross to QTL results", "years"),
    c("generations", "minimize", "Crossing and selfing generations", "generations")
  )
  out <- lapply(rows, function(r) new_sim_metric(r[1], r[2], label = r[3], units = r[4]))
  stats::setNames(out, vapply(rows, `[`, character(1), 1))
}

new_sim_metric <- function(name, direction, fn = NULL, label = name, units = "") {
  direction <- rlang::arg_match0(direction, c("maximize", "minimize", "zero"))
  structure(list(name = name, direction = direction, fn = fn, label = label, units = units),
            class = "sim_metric")
}

#' Define a simulation metric
#'
#' A metric reads one number from a replicate's outcomes and says whether
#' bigger or smaller is better. Built-in metrics (see [sim_metric_set()]) need
#' no definition; use `sim_metric()` for anything else.
#'
#' @param name Metric name; also the outcome to read when `fn` is `NULL`.
#' @param direction `"maximize"`, `"minimize"` or `"zero"` (closest to zero is best).
#' @param fn Optional function of the outcome list returning one number, for a derived metric.
#' @param label,units Labels used in printing and plots.
#' @return A `sim_metric` object.
#' @examples
#' cost_per_line <- sim_metric("cost_per_line", "minimize", fn = function(o) o$cost / o$lines)
#' sim_metric_set(power, cost, cost_per_line)
#' @export
sim_metric <- function(name, direction = "maximize", fn = NULL, label = name, units = "") {
  new_sim_metric(name, direction, fn, label, units)
}

#' Combine metrics into a metric set
#'
#' Like `yardstick::metric_set()`. Name built-in metrics bare or as strings and
#' add custom ones made with [sim_metric()]. The first metric is the one tuners
#' optimize and selection helpers rank by unless told otherwise.
#'
#' Built-in metrics: `power`, `false_positive`, `location_error`,
#' `support_interval`, `peak_lod`, `realized_pve`, `threshold` (mapping);
#' `release_gain`, `gain_per_year`, `gain_per_cost`, `expected_heterozygosity`,
#' `genetic_variance` (breeding programs); `lines`, `plots`, `genotyped`,
#' `crosses`, `cost`, `years`, `generations` (resources).
#'
#' @param ... Metric names (bare or quoted) or `sim_metric` objects.
#' @return A `sim_metric_set`.
#' @examples
#' sim_metric_set(power, cost, years)
#' @export
sim_metric_set <- function(...) {
  quos <- rlang::enquos(...)
  if (!length(quos)) abort("Give at least one metric.")
  known <- builtin_metrics()
  metrics <- lapply(quos, function(q) {
    expr <- rlang::quo_get_expr(q)
    if (rlang::is_symbol(expr) && rlang::as_string(expr) %in% names(known)) {
      if (!exists(rlang::as_string(expr), envir = rlang::quo_get_env(q)) ||
          !inherits(rlang::eval_tidy(q), "sim_metric")) {
        return(known[[rlang::as_string(expr)]])
      }
    }
    value <- if (rlang::is_symbol(expr)) {
      tryCatch(rlang::eval_tidy(q), error = function(e) rlang::as_string(expr))
    } else {
      rlang::eval_tidy(q)
    }
    as_sim_metric(value, known)
  })
  names(metrics) <- vapply(metrics, function(m) m$name, character(1))
  if (anyDuplicated(names(metrics))) abort("Metric names must be unique.")
  structure(metrics, class = "sim_metric_set")
}

as_sim_metric <- function(x, known = builtin_metrics()) {
  if (inherits(x, "sim_metric")) return(x)
  if (is.character(x) && length(x) == 1) {
    if (!x %in% names(known)) {
      abort(c(sprintf("Unknown metric `%s`.", x), i = "Define it with sim_metric()."))
    }
    return(known[[x]])
  }
  abort("Metrics must be names or sim_metric() objects.")
}

as_metric_set <- function(metrics, spec) {
  if (is.null(metrics)) {
    metrics <- sim_default_metrics(spec)
    if (is.null(metrics)) abort("Give a metric set with sim_metric_set().")
  }
  if (inherits(metrics, "sim_metric_set")) return(metrics)
  if (inherits(metrics, "sim_metric")) metrics <- list(metrics)
  out <- lapply(metrics, as_sim_metric)
  names(out) <- vapply(out, function(m) m$name, character(1))
  structure(out, class = "sim_metric_set")
}

metric_value <- function(m, outcome) {
  if (!is.null(m$fn)) return(as.numeric(m$fn(outcome)))
  v <- outcome[[m$name]]
  if (is.null(v)) {
    abort(sprintf("The simulation did not report `%s`; it reported %s.", m$name,
                  paste(names(outcome), collapse = ", ")))
  }
  as.numeric(v)
}

metric_sign <- function(m) if (identical(m$direction, "maximize")) 1 else -1

#' @export
print.sim_metric_set <- function(x, ...) {
  cat("A simulation metric set:\n")
  for (m in x) cat(sprintf("  %-18s %-9s %s\n", m$name, m$direction, m$label))
  invisible(x)
}

#' @export
print.sim_metric <- function(x, ...) {
  cat(sprintf("Simulation metric %s (%s): %s\n", x$name, x$direction, x$label))
  invisible(x)
}
