# Tidy results and selection helpers -----------------------------------------------------------

summarize_raw <- function(raw, metrics, by_scenario = FALSE) {
  if (is.null(raw) || !nrow(raw)) {
    return(tibble::tibble(.config = character(), .metric = character(), mean = numeric(), n = integer(),
                          std_err = numeric()))
  }
  groups <- if (by_scenario) paste(raw$.config, raw$.scenario, sep = "\r") else raw$.config
  out <- list()
  for (g in unique(groups)) {
    sub <- raw[groups == g, , drop = FALSE]
    for (m in metrics) {
      # Average scenarios within each seed first, so standard errors stay paired by seed.
      v <- if (by_scenario) sub[[m$name]] else as.numeric(tapply(sub[[m$name]], sub$seed, mean))
      v <- v[is.finite(v)]
      n <- length(v)
      row <- tibble::tibble(.config = sub$.config[1], .metric = m$name, direction = m$direction,
                            mean = if (n) mean(v) else NA_real_, n = n,
                            std_err = if (n > 1) stats::sd(v) / sqrt(n) else NA_real_)
      if (by_scenario) row$.scenario <- sub$.scenario[1]
      out[[length(out) + 1]] <- row
    }
  }
  do.call(rbind, out)
}

constraint_rule <- function(rule, metric) {
  if (is.character(rule) && length(rule) == 2) return(list(op = rule[1], bound = as.numeric(rule[2])))
  if (is.list(rule) && length(rule) == 2) return(list(op = rule[[1]], bound = as.numeric(rule[[2]])))
  if (!is.numeric(rule) || length(rule) != 1) abort("Each constraint is a number or c(op, value).")
  list(op = if (identical(metric$direction, "maximize")) ">=" else "<=", bound = rule)
}

meets_constraints <- function(summ, constraints, metrics) {
  for (name in names(constraints)) {
    if (!name %in% names(metrics)) abort(sprintf("Constraint `%s` is not in the metric set.", name))
    rule <- constraint_rule(constraints[[name]], metrics[[name]])
    v <- summ$mean[summ$.metric == name]
    if (!length(v) || !is.finite(v)) return(FALSE)
    ok <- switch(rule$op, ">=" = v >= rule$bound, ">" = v > rule$bound, "<=" = v <= rule$bound,
                 "<" = v < rule$bound, abort(sprintf("Unknown constraint operator %s.", rule$op)))
    if (!ok) return(FALSE)
  }
  TRUE
}

best_configs <- function(summ, metric, constraints, metrics, n = 5) {
  configs <- unique(summ$.config)
  if (length(constraints)) {
    configs <- configs[vapply(configs, function(cfg) meets_constraints(summ[summ$.config == cfg, ], constraints,
                                                                       metrics), logical(1))]
  }
  s <- summ[summ$.metric == metric$name & summ$.config %in% configs & is.finite(summ$mean), , drop = FALSE]
  key <- if (identical(metric$direction, "zero")) abs(s$mean) else -metric_sign(metric) * s$mean
  utils::head(s[order(key), , drop = FALSE], n)
}

check_results <- function(x) {
  if (!inherits(x, "sim_results")) abort("`x` must be results from a tidyagsim tuner.")
}

get_metric <- function(x, metric) {
  if (is.null(metric)) return(x$metrics[[1]])
  if (!metric %in% names(x$metrics)) abort(sprintf("`%s` is not in the metric set.", metric))
  x$metrics[[metric]]
}

# collect_metrics ---------------------------------------------------------------------------------

#' Collect tidy metrics
#'
#' Like `tune::collect_metrics()`. With `summarize = TRUE`, one row per design
#' (and scenario, when there are several) and metric, with the mean across
#' replicates, the number of replicates `n` and the Monte Carlo standard error
#' `std_err`. When scenarios are pooled, each seed's value is first averaged
#' over scenarios so the standard error stays paired by seed.
#'
#' @param x A `sim_results` object.
#' @param ... Not used.
#' @param summarize Summarize over replicates (`TRUE`) or return one row per
#'   design, scenario and replicate (`FALSE`).
#' @param by_scenario Summarize within each scenario. Defaults to `TRUE` when
#'   the study has more than one scenario.
#' @param type `"long"` (one row per metric, like tune) or `"wide"` (one row per
#'   design, a column per metric and `<metric>_se` columns).
#' @return A tibble.
#' @export
collect_metrics <- function(x, ...) UseMethod("collect_metrics")

#' @export
collect_metrics.default <- function(x, ...) forward("tune", "collect_metrics", x, ...)

#' @rdname collect_metrics
#' @export
collect_metrics.sim_results <- function(x, ..., summarize = TRUE, by_scenario = NULL,
                                        type = c("long", "wide")) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  ids <- x$params$id
  if (!summarize) {
    return(x$raw[, setdiff(names(x$raw), ".outcome"), drop = FALSE])
  }
  by_scenario <- by_scenario %||% (nrow(attr(x$resamples, "scenarios")) > 1)
  summ <- summarize_raw(x$raw, x$metrics, by_scenario)
  values <- x$configs[match(summ$.config, x$configs$.config), ids, drop = FALSE]
  long <- tibble::as_tibble(cbind(values, summ[, setdiff(names(summ), ".config"), drop = FALSE],
                                  summ[, ".config", drop = FALSE]))
  if (by_scenario) long <- long[, c(setdiff(names(long), c(".scenario", ".config")), ".scenario", ".config")]
  if (type == "long") return(long)
  key <- if (by_scenario) c(".config", ".scenario") else ".config"
  designs <- unique(long[, c(ids, key), drop = FALSE])
  for (m in names(x$metrics)) {
    sub <- long[long$.metric == m, , drop = FALSE]
    k1 <- do.call(paste, c(designs[key], sep = "\r"))
    k2 <- do.call(paste, c(sub[key], sep = "\r"))
    designs[[m]] <- sub$mean[match(k1, k2)]
    designs[[paste0(m, "_se")]] <- sub$std_err[match(k1, k2)]
  }
  designs$n <- long$n[long$.metric == names(x$metrics)[1]][match(
    do.call(paste, c(designs[key], sep = "\r")),
    do.call(paste, c(long[long$.metric == names(x$metrics)[1], key], sep = "\r")))]
  designs[, c(ids, setdiff(names(designs), c(ids, key)), key)]
}

# show_best / select_* ----------------------------------------------------------------------------

#' Find the best designs
#'
#' `show_best()` lists the top designs by the mean of one metric;
#' `select_best()` returns the parameter values of the best one;
#' `select_by_one_std_err()` returns the simplest design within one standard
#' error of the best, and `select_by_pct_loss()` the simplest within a
#' percentage loss. These mirror the tune functions of the same names.
#'
#' Means are pooled over scenarios unless `scenario` names one. Constraints
#' restrict which designs are eligible, e.g. `constraints = list(power = 0.8)`
#' (at least 0.8, because power is maximized) or
#' `list(cost = c("<=", 20000))`.
#'
#' @param x A `sim_results` object.
#' @param ... For `select_by_one_std_err()` and `select_by_pct_loss()`: columns
#'   (parameters or metrics) that order designs from simplest, e.g. `lines` or
#'   `desc(families)` or `cost`.
#' @param metric Metric to rank by; defaults to the first in the set.
#' @param n Number of designs to show.
#' @param constraints Named list of bounds on metric means.
#' @param scenario A scenario label to rank within, or `NULL` to pool.
#' @param limit Allowed percentage loss for `select_by_pct_loss()`.
#' @return A tibble.
#' @export
show_best <- function(x, ...) UseMethod("show_best")

#' @export
show_best.default <- function(x, ...) forward("tune", "show_best", x, ...)

#' @rdname show_best
#' @export
show_best.sim_results <- function(x, ..., metric = NULL, n = 5, constraints = NULL, scenario = NULL) {
  rlang::check_dots_empty()
  m <- get_metric(x, metric)
  summ <- scenario_summary(x, scenario)
  best <- best_configs(summ, m, constraints, x$metrics, n)
  values <- x$configs[match(best$.config, x$configs$.config), x$params$id, drop = FALSE]
  tibble::as_tibble(cbind(values, best[, c(".metric", "mean", "n", "std_err")], best[, ".config", drop = FALSE]))
}

scenario_summary <- function(x, scenario) {
  if (is.null(scenario)) return(summarize_raw(x$raw, x$metrics, by_scenario = FALSE))
  if (!scenario %in% x$raw$.scenario) abort(sprintf("No scenario labelled `%s`.", scenario))
  summarize_raw(x$raw[x$raw$.scenario == scenario, , drop = FALSE], x$metrics, by_scenario = FALSE)
}

#' @rdname show_best
#' @export
select_best <- function(x, ...) UseMethod("select_best")

#' @export
select_best.default <- function(x, ...) forward("tune", "select_best", x, ...)

#' @rdname show_best
#' @export
select_best.sim_results <- function(x, ..., metric = NULL, constraints = NULL, scenario = NULL) {
  rlang::check_dots_empty()
  best <- show_best.sim_results(x, metric = metric, n = 1, constraints = constraints, scenario = scenario)
  if (!nrow(best)) abort("No design meets the constraints.")
  best[, c(x$params$id, ".config")]
}

#' @rdname show_best
#' @export
select_by_one_std_err <- function(x, ...) UseMethod("select_by_one_std_err")

#' @export
select_by_one_std_err.default <- function(x, ...) forward("tune", "select_by_one_std_err", x, ...)

#' @rdname show_best
#' @export
select_by_one_std_err.sim_results <- function(x, ..., metric = NULL, constraints = NULL, scenario = NULL) {
  select_within(x, rlang::enquos(...), metric, constraints, scenario, "se", 1)
}

#' @rdname show_best
#' @export
select_by_pct_loss <- function(x, ...) UseMethod("select_by_pct_loss")

#' @export
select_by_pct_loss.default <- function(x, ...) forward("tune", "select_by_pct_loss", x, ...)

#' @rdname show_best
#' @export
select_by_pct_loss.sim_results <- function(x, ..., metric = NULL, limit = 2, constraints = NULL, scenario = NULL) {
  select_within(x, rlang::enquos(...), metric, constraints, scenario, "pct", limit)
}

select_within <- function(x, order_by, metric, constraints, scenario, rule, amount) {
  if (!length(order_by)) abort("Name at least one column that orders designs from simplest, e.g. `lines`.")
  m <- get_metric(x, metric)
  summ <- scenario_summary(x, scenario)
  best <- best_configs(summ, m, constraints, x$metrics, 1)
  if (!nrow(best)) abort("No design meets the constraints.")
  eligible <- best_configs(summ, m, constraints, x$metrics, Inf)
  bound <- if (rule == "se") {
    se <- if (is.finite(best$std_err)) best$std_err else 0
    best$mean - metric_sign(m) * se
  } else {
    best$mean - metric_sign(m) * abs(best$mean) * amount / 100
  }
  near <- eligible[metric_sign(m) * (eligible$mean - bound) >= -1e-12, , drop = FALSE]
  table <- x$configs[match(near$.config, x$configs$.config), , drop = FALSE]
  wide <- collect_metrics.sim_results(x, by_scenario = FALSE, type = "wide")
  for (name in setdiff(names(x$metrics), names(table))) table[[name]] <- wide[[name]][match(table$.config, wide$.config)]
  keys <- lapply(order_by, function(q) rlang::eval_tidy(q, data = table, env = rlang::env(desc = function(v) {
    if (is.character(v) || is.factor(v)) -xtfrm(v) else -v
  })))
  table <- table[do.call(order, keys), , drop = FALSE]
  table[1, c(x$params$id, ".config")]
}

# Pareto fronts and paired comparisons -------------------------------------------------------------

#' Pareto front of designs
#'
#' Designs that no other design beats on every listed metric (comparing means,
#' respecting each metric's direction). The default front is the first two
#' metrics, for example power against cost.
#'
#' @param x A `sim_results` object.
#' @param metrics Metric names; defaults to the first two in the set.
#' @param constraints,scenario As in [show_best()].
#' @return A tibble with one row per design on the front, with means and standard errors.
#' @export
pareto_front <- function(x, metrics = NULL, constraints = NULL, scenario = NULL) {
  check_results(x)
  metrics <- metrics %||% utils::head(names(x$metrics), 2)
  ms <- lapply(metrics, get_metric, x = x)
  summ <- scenario_summary(x, scenario)
  configs <- unique(summ$.config)
  if (length(constraints)) {
    configs <- configs[vapply(configs, function(cfg) meets_constraints(summ[summ$.config == cfg, ], constraints,
                                                                       x$metrics), logical(1))]
  }
  pts <- sapply(ms, function(m) {
    s <- summ[summ$.metric == m$name, ]
    metric_sign(m) * s$mean[match(configs, s$.config)]
  })
  pts <- matrix(pts, nrow = length(configs))
  keep <- vapply(seq_along(configs), function(i) {
    if (any(!is.finite(pts[i, ]))) return(FALSE)
    !any(vapply(seq_along(configs), function(j) {
      j != i && all(is.finite(pts[j, ])) && all(pts[j, ] >= pts[i, ]) && any(pts[j, ] > pts[i, ])
    }, logical(1)))
  }, logical(1))
  front <- configs[keep]
  wide <- collect_metrics.sim_results(x, by_scenario = FALSE, type = "wide")
  if (!is.null(scenario)) {
    wide <- collect_metrics.sim_results(x, by_scenario = TRUE, type = "wide")
    wide <- wide[wide$.scenario == scenario, , drop = FALSE]
  }
  out <- wide[match(front, wide$.config), c(x$params$id, as.vector(rbind(metrics, paste0(metrics, "_se"))),
                                            ".config"), drop = FALSE]
  out[order(out[[metrics[1]]]), , drop = FALSE]
}

#' Paired comparison of two designs
#'
#' The difference `b - a` in a metric on shared seeds (and scenarios), with
#' its standard error and the one-sided p-value that `b` is worse than `a`.
#' Pairing on seeds is what makes small differences between designs detectable.
#'
#' @param x A `sim_results` object.
#' @param a,b `.config` labels.
#' @param metric Metric name; defaults to the first in the set.
#' @return A one-row tibble.
#' @export
compare_designs <- function(x, a, b, metric = NULL) {
  check_results(x)
  m <- get_metric(x, metric)
  ra <- x$raw[x$raw$.config == a, ]
  rb <- x$raw[x$raw$.config == b, ]
  ka <- paste(ra$.scenario, ra$seed)
  kb <- paste(rb$.scenario, rb$seed)
  shared <- intersect(ka, kb)
  d <- rb[[m$name]][match(shared, kb)] - ra[[m$name]][match(shared, ka)]
  d <- d[is.finite(d)]
  tibble::tibble(.metric = m$name, a = a, b = b, mean_difference = mean(d), std_err = stats::sd(d) / sqrt(length(d)),
                 n = length(d), p_b_worse = paired_worse_p(metric_sign(m) * d))
}

#' @export
print.sim_results <- function(x, ...) {
  sc <- attr(x$resamples, "scenarios")
  cat(sprintf("# Simulation study (%s): %d designs x %d replicates x %d scenario(s) = %d simulations\n",
              x$method, nrow(x$configs), nrow(x$resamples), nrow(sc), nrow(x$raw)))
  cat(sprintf("# Spec: %s; metrics: %s\n", class(x$spec)[1], paste(names(x$metrics), collapse = ", ")))
  if (nrow(x$infeasible)) cat(sprintf("# %d design/scenario pair(s) infeasible and skipped\n", nrow(x$infeasible)))
  if (nrow(x$notes)) cat(sprintf("# %d batch(es) failed; see x$notes\n", nrow(x$notes)))
  if (nrow(x$raw)) {
    cat(sprintf("# Best by %s:\n", names(x$metrics)[1]))
    print(show_best.sim_results(x, n = 3))
  }
  invisible(x)
}

#' Notes from failed simulations
#' @param x A `sim_results` object.
#' @return A tibble of failures and a tibble of infeasible designs, in a list.
#' @export
collect_notes <- function(x) {
  check_results(x)
  list(failures = x$notes, infeasible = x$infeasible)
}
