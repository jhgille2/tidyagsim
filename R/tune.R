# Tuners -------------------------------------------------------------------------------------

with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }
  })
  # Fix the generator too: parallel workers (future.seed = TRUE) switch to
  # L'Ecuyer-CMRG, which would give different draws from the same seed.
  # Restoring .Random.seed afterwards restores the caller's generator.
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}

param_ids <- function(params) params$id

param_size <- function(p) {
  if (inherits(p, "qual_param")) return(length(p$values))
  if (identical(p$type, "integer") && is.null(p$trans)) return(diff(unlist(p$range)) + 1)
  Inf
}

space_filling <- function(params, size) {
  ns <- asNamespace("dials")
  if (exists("grid_space_filling", envir = ns)) {
    return(getExportedValue("dials", "grid_space_filling")(params, size = size, type = "max_entropy"))
  }
  if (all(vapply(params$object, inherits, logical(1), "quant_param"))) {
    return(dials::grid_max_entropy(params, size = size))
  }
  dials::grid_latin_hypercube(params, size = size)
}

#' Build the candidate grid for a study
#'
#' `NULL` or a number gives that many space-filling designs (10 by default),
#' or every combination when the space is that small. A data frame is used as
#' given after checking its columns. `expand_grid`-style tibbles from
#' [dials::grid_regular()] work directly.
#'
#' @param params A dials `parameters` object.
#' @param grid `NULL`, a number of designs, or a data frame of designs.
#' @param seed Seed for the space-filling design.
#' @return A tibble with `.config` and one column per parameter id.
#' @keywords internal
#' @export
sim_grid <- function(params, grid = NULL, seed = 1) {
  ids <- param_ids(params)
  if (is.null(grid) || (is.numeric(grid) && length(grid) == 1)) {
    size <- grid %||% 10
    total <- prod(vapply(params$object, param_size, numeric(1)))
    grid <- with_seed(seed, {
      if (is.finite(total) && total <= size) {
        levels <- vapply(params$object, function(p) as.integer(param_size(p)), integer(1))
        dials::grid_regular(params, levels = levels)
      } else {
        space_filling(params, size)
      }
    })
    grid <- grid[do.call(order, unname(as.list(grid))), , drop = FALSE]
  }
  grid <- tibble::as_tibble(grid)
  missing <- setdiff(ids, names(grid))
  extra <- setdiff(names(grid), c(ids, ".config"))
  if (length(missing) || length(extra)) {
    abort(sprintf("Grid columns must match the tuned parameter ids (%s). Missing: %s. Extra: %s.",
                  paste(ids, collapse = ", "), paste(missing, collapse = ", ") %|e|% "none",
                  paste(extra, collapse = ", ") %|e|% "none"))
  }
  grid <- unique(grid[, ids, drop = FALSE])
  width <- max(3, nchar(nrow(grid)))
  tibble::add_column(grid, .config = sprintf(paste0("Config%0", width, "d"), seq_len(nrow(grid))), .before = 1)
}

`%|e|%` <- function(a, b) if (!length(a) || identical(a, "")) b else a

prepare_study <- function(object, resamples, metrics, param_info, grid = NULL) {
  if (!inherits(object, "sim_spec")) abort("`object` must be a simulation spec.")
  resamples <- as_rset(resamples)
  metrics <- as_metric_set(metrics, object)
  params <- param_info %||% extract_parameter_set_dials.sim_spec(object, grid = grid)
  marked <- tune_args.sim_spec(object)
  if (!setequal(params$id, unique(marked$id))) {
    abort(sprintf("`param_info` ids (%s) must match the tuned parameters (%s).",
                  paste(params$id, collapse = ", "), paste(unique(marked$id), collapse = ", ")))
  }
  list(resamples = resamples, metrics = metrics, params = params)
}

new_sim_results <- function(object, study, configs, rows, runner, method, trace = NULL) {
  raw <- bind_rows_list(rows)
  ids <- study$params$id
  if (is.null(raw)) {
    raw <- tibble::tibble(.config = character(), .scenario = character(), id = character(), seed = integer())
    for (m in study$metrics) raw[[m$name]] <- numeric()
  }
  raw <- merge_configs(raw, configs, ids)
  ord <- order(match(raw$.config, configs$.config), raw$.scenario, raw$seed)
  raw <- raw[ord, , drop = FALSE]
  notes <- if (length(runner$notes)) do.call(rbind, runner$notes) else
    tibble::tibble(.config = character(), .scenario = character(), note = character())
  infeasible <- if (length(runner$infeasible)) do.call(rbind, unname(runner$infeasible)) else
    tibble::tibble(.config = character(), .scenario = character(), reason = character())
  if (nrow(notes)) {
    warn(sprintf("%d design/scenario batch(es) failed. See `x$notes`; first: %s", nrow(notes), notes$note[1]))
  }
  structure(list(raw = raw, configs = configs, params = study$params, metrics = study$metrics,
                 resamples = study$resamples, spec = object, method = method, notes = notes,
                 infeasible = infeasible, trace = trace %||% tibble::tibble()),
            class = c(paste0("sim_results_", method), "sim_results"))
}

merge_configs <- function(raw, configs, ids) {
  idx <- match(raw$.config, configs$.config)
  values <- configs[idx, ids, drop = FALSE]
  tibble::as_tibble(cbind(raw[, ".config", drop = FALSE], values, raw[, setdiff(names(raw), ".config"), drop = FALSE]))
}

# Generics that forward to tune when given tidymodels objects ---------------------------------

#' Run every design in a grid
#'
#' The simulation analogue of `tune::tune_grid()`: every design runs on every
#' replicate seed and scenario.
#'
#' @param object A simulation spec with `tune()` arguments.
#' @param resamples Replicates from [sim_replicates()], or a number of replicates.
#' @param ... Not used.
#' @param grid `NULL` (10 space-filling designs), a number of designs, or a
#'   data frame of designs with one column per tuned parameter id.
#' @param metrics A [sim_metric_set()]; defaults to the spec's own.
#' @param param_info A dials `parameters` object with changed ranges, from
#'   `extract_parameter_set_dials(spec) |> update(...)`.
#' @param control A [control_sim()] object.
#' @param seed Seed for generating the space-filling grid.
#' @return A `sim_results` object; see [collect_metrics()].
#' @examples
#' spec <- mapping_population(type = "RIL", lines = tune(), n_chromosomes = 2, n_perm = 20)
#' res <- tune_grid(spec, sim_replicates(3), grid = data.frame(lines = c(60, 120)),
#'                  metrics = sim_metric_set(power, cost))
#' collect_metrics(res)
#' @export
tune_grid <- function(object, ...) UseMethod("tune_grid")

#' @export
tune_grid.default <- function(object, ...) forward("tune", "tune_grid", object, ...)

#' @rdname tune_grid
#' @export
tune_grid.sim_spec <- function(object, resamples = NULL, ..., grid = NULL, metrics = NULL, param_info = NULL,
                               control = control_sim(), seed = 1) {
  rlang::check_dots_empty()
  study <- prepare_study(object, resamples, metrics, param_info, grid)
  if (!nrow(study$params)) abort("Mark at least one argument with tune(), or use fit_resamples().")
  configs <- sim_grid(study$params, grid, seed)
  runner <- new_runner(object, study$metrics, control)
  rows <- run_designs(runner, configs, scenario_list(study$resamples), study$resamples)
  new_sim_results(object, study, configs, rows, runner, "grid")
}

#' Run one fixed design on every replicate
#'
#' The simulation analogue of `tune::fit_resamples()`.
#'
#' @inheritParams tune_grid
#' @return A `sim_results` object with one design.
#' @examples
#' spec <- mapping_population(type = "F2", lines = 100, n_chromosomes = 2, n_perm = 20)
#' collect_metrics(fit_resamples(spec, 3, metrics = sim_metric_set(power)))
#' @export
fit_resamples <- function(object, ...) UseMethod("fit_resamples")

#' @export
fit_resamples.default <- function(object, ...) forward("tune", "fit_resamples", object, ...)

#' @rdname fit_resamples
#' @export
fit_resamples.sim_spec <- function(object, resamples = NULL, ..., metrics = NULL, control = control_sim()) {
  rlang::check_dots_empty()
  if (nrow(tune_args.sim_spec(object))) abort("Fix or finalize the tuned arguments first (see finalize_sim()).")
  study <- prepare_study(object, resamples, metrics, NULL)
  configs <- tibble::tibble(.config = "Config001")
  runner <- new_runner(object, study$metrics, control)
  rows <- run_designs(runner, configs, scenario_list(study$resamples), study$resamples)
  new_sim_results(object, study, configs, rows, runner, "resamples")
}

# Racing ----------------------------------------------------------------------------------------

paired_worse_p <- function(diffs) {
  n <- length(diffs)
  if (n < 2) return(1)
  s <- stats::sd(diffs)
  m <- mean(diffs)
  if (!is.finite(s) || s < 1e-12) return(if (m < 0) 0 else 1)
  stats::pt(m / (s / sqrt(n)), df = n - 1)
}

seed_scores <- function(raw, metric) {
  # Mean over scenarios for each config and seed, oriented so bigger is better.
  key <- paste(raw$.config, raw$seed, sep = "\r")
  v <- tapply(metric_sign(metric) * raw[[metric$name]], key, mean)
  parts <- do.call(rbind, strsplit(names(v), "\r", fixed = TRUE))
  tibble::tibble(.config = parts[, 1], seed = as.integer(parts[, 2]), score = as.numeric(v))
}

#' Race designs and drop the clearly worse ones early
#'
#' Like `finetune::tune_race_anova()`. Every design runs on the first `burn_in`
#' replicates. After each further replicate, a design is dropped when a
#' one-sided paired t-test on shared seeds says it is worse than the current
#' leader at level `alpha`. Designs that break `constraints` on current means
#' are dropped as soon as a design meeting them leads. Survivors are run on
#' every replicate, so their estimates are comparable with a full grid.
#'
#' @inheritParams tune_grid
#' @param burn_in Replicates every design gets before any is dropped (at least 2).
#' @param alpha Significance level for dropping a design.
#' @param metric The metric to race on; defaults to the first in the set.
#' @param constraints Named list of bounds on other metrics' means, e.g.
#'   `list(power = 0.8)`. A bare number is a floor for metrics to maximize and
#'   a ceiling for metrics to minimize; `list(cost = c("<=", 5000))` is explicit.
#' @return A `sim_results` object; `x$trace` records how many designs survived each round.
#' @export
tune_race <- function(object, resamples = NULL, ..., grid = NULL, metrics = NULL, param_info = NULL,
                      control = control_sim(), burn_in = 3, alpha = 0.05, metric = NULL, constraints = NULL,
                      seed = 1) {
  rlang::check_dots_empty()
  study <- prepare_study(object, resamples, metrics, param_info, grid)
  reps <- study$resamples
  if (burn_in < 2 || burn_in > nrow(reps)) abort("`burn_in` must be between 2 and the number of replicates.")
  m <- if (is.null(metric)) study$metrics[[1]] else study$metrics[[metric]]
  configs <- sim_grid(study$params, grid, seed)
  scen <- scenario_list(reps)
  runner <- new_runner(object, study$metrics, control)
  rows <- run_designs(runner, configs, scen, reps[seq_len(burn_in), ])
  alive <- configs
  trace <- list()
  for (step in burn_in:nrow(reps)) {
    raw <- bind_rows_list(rows)
    if (is.null(raw)) break
    alive <- alive[alive$.config %in% raw$.config, , drop = FALSE]
    if (nrow(alive) <= 1) break
    sc <- seed_scores(raw[raw$.config %in% alive$.config, ], m)
    means <- tapply(sc$score, sc$.config, mean)
    ok <- rep(TRUE, length(means))
    names(ok) <- names(means)
    if (length(constraints)) {
      summ <- summarize_raw(raw[raw$.config %in% alive$.config, ], study$metrics, by_scenario = FALSE)
      ok <- vapply(names(means), function(cfg) meets_constraints(summ[summ$.config == cfg, ], constraints,
                                                                 study$metrics), logical(1))
    }
    ranked <- names(means)[order(!ok, -means)]
    leader <- ranked[1]
    keep <- vapply(alive$.config, function(cfg) {
      if (cfg == leader) return(TRUE)
      if (length(constraints) && ok[[leader]] && !ok[[cfg]]) return(FALSE)
      a <- sc[sc$.config == cfg, ]
      b <- sc[sc$.config == leader, ]
      shared <- intersect(a$seed, b$seed)
      diffs <- a$score[match(shared, a$seed)] - b$score[match(shared, b$seed)]
      paired_worse_p(diffs) >= alpha
    }, logical(1))
    trace[[length(trace) + 1]] <- tibble::tibble(replicates = step, designs = nrow(alive), kept = sum(keep),
                                                 leader = leader)
    if (isTRUE(control$verbose)) {
      cli::cli_inform("After {step} replicate{?s}: {sum(keep)} of {nrow(alive)} designs remain (leader {leader}).")
    }
    alive <- alive[keep, , drop = FALSE]
    if (step < nrow(reps) && nrow(alive) > 1) {
      rows <- c(rows, run_designs(runner, alive, scen, reps[step + 1, ]))
    }
  }
  # Finish survivors on the whole panel; cached replicates are not rerun.
  done <- bind_rows_list(rows)
  finish <- run_designs(runner, alive, scen, reps)
  have <- paste(done$.config, done$.scenario, done$seed)
  rows <- c(rows, Filter(function(r) !paste(r$.config, r$.scenario, r$seed) %in% have, finish))
  new_sim_results(object, study, configs, rows, runner, "race",
                  if (length(trace)) do.call(rbind, trace) else NULL)
}

# Bayesian optimization ---------------------------------------------------------------------------

encode_values <- function(params, values) {
  cols <- lapply(seq_len(nrow(params)), function(i) {
    p <- params$object[[i]]
    v <- values[[params$id[i]]]
    if (inherits(p, "qual_param")) {
      return(vapply(p$values, function(level) as.numeric(as.character(v) == as.character(level)), numeric(length(v))))
    }
    rng <- unlist(p$range)
    x <- as.numeric(v)
    if (!is.null(p$trans)) {
      x <- p$trans$transform(x)
    }
    if (diff(rng) > 0) (x - rng[1]) / diff(rng) else x * 0
  })
  do.call(cbind, lapply(cols, as.matrix))
}

gp_fit <- function(x, y, noise) {
  mu0 <- mean(y)
  s0 <- max(stats::sd(y), 1e-9)
  if (!is.finite(s0)) s0 <- 1
  z <- (y - mu0) / s0
  nz <- noise / s0^2
  d <- as.matrix(stats::dist(x))
  matern <- function(d, ell) {
    r <- sqrt(5) * d / ell
    (1 + r + r^2 / 3) * exp(-r)
  }
  best <- NULL
  for (ell in c(0.1, 0.2, 0.35, 0.5, 0.8, 1.2, 2)) {
    k <- matern(d, ell) + diag(nz + 1e-6, length(z))
    ch <- tryCatch(chol(k), error = function(e) NULL)
    if (is.null(ch)) next
    alpha <- backsolve(ch, forwardsolve(t(ch), z))
    lml <- -0.5 * sum(z * alpha) - sum(log(diag(ch)))
    if (is.null(best) || lml > best$lml) best <- list(lml = lml, ell = ell, ch = ch, alpha = alpha)
  }
  function(xs) {
    ds <- sqrt(pmax(outer(rowSums(xs^2), rowSums(x^2), "+") - 2 * xs %*% t(x), 0))
    ks <- matern(ds, best$ell)
    mu <- as.vector(ks %*% best$alpha)
    v <- forwardsolve(t(best$ch), t(ks))
    sdv <- sqrt(pmax(1 - colSums(v^2), 1e-12))
    list(mean = mu * s0 + mu0, sd = sdv * s0)
  }
}

expected_improvement <- function(mu, sd, best) {
  z <- (mu - best) / sd
  (mu - best) * stats::pnorm(z) + sd * stats::dnorm(z)
}

#' Bayesian optimization of a simulation study
#'
#' Like `tune::tune_bayes()`. Starts from `initial` space-filling designs, then
#' proposes one design per iteration by expected improvement under a Gaussian
#' process fitted to design means (with Monte Carlo standard errors as noise).
#' With `constraints`, improvement is weighted by the modeled probability that
#' each constrained metric meets its bound, so "the cheapest design with at
#' least 80% power" is answered directly: `metrics = sim_metric_set(cost, power)`,
#' `constraints = list(power = 0.8)`.
#'
#' @inheritParams tune_race
#' @param initial Number of starting designs, or a data frame of them.
#' @param iter Number of proposals after the initial designs.
#' @param no_improve Stop after this many iterations without a new best design.
#' @param candidates Size of the random candidate pool scored at each iteration.
#' @return A `sim_results` object; `x$trace` lists each proposal.
#' @export
tune_bayes <- function(object, ...) UseMethod("tune_bayes")

#' @export
tune_bayes.default <- function(object, ...) forward("tune", "tune_bayes", object, ...)

#' @rdname tune_bayes
#' @export
tune_bayes.sim_spec <- function(object, resamples = NULL, ..., metrics = NULL, param_info = NULL, initial = 6,
                                iter = 10, metric = NULL, constraints = NULL, no_improve = Inf,
                                candidates = 2000, control = control_sim(), seed = 1) {
  rlang::check_dots_empty()
  study <- prepare_study(object, resamples, metrics, param_info)
  if (!nrow(study$params)) abort("Mark at least one argument with tune().")
  m <- if (is.null(metric)) study$metrics[[1]] else study$metrics[[metric]]
  configs <- sim_grid(study$params, initial, seed)
  scen <- scenario_list(study$resamples)
  runner <- new_runner(object, study$metrics, control)
  rows <- run_designs(runner, configs, scen, study$resamples)
  ids <- study$params$id
  pool <- with_seed(seed + 1, dials::grid_random(study$params, size = candidates))
  pool <- unique(pool)
  key <- function(df) do.call(paste, c(lapply(ids, function(i) format(df[[i]], digits = 10)), sep = "\r"))
  trace <- list()
  best_so_far <- NA_character_
  stale <- 0
  for (it in seq_len(iter)) {
    raw <- merge_configs(bind_rows_list(rows), configs, ids)
    summ <- summarize_raw(raw, study$metrics, by_scenario = FALSE)
    observed <- configs[configs$.config %in% summ$.config, , drop = FALSE]
    if (nrow(observed) < 2) break
    stat_of <- function(name, what) {
      s <- summ[summ$.metric == name, ]
      s[[what]][match(observed$.config, s$.config)]
    }
    x <- encode_values(study$params, observed)
    model <- function(metric) {
      y <- metric_sign(metric) * stat_of(metric$name, "mean")
      se <- stat_of(metric$name, "std_err")
      se[!is.finite(se)] <- if (any(is.finite(se))) mean(se[is.finite(se)]) else 0
      ok <- is.finite(y)
      list(fit = gp_fit(x[ok, , drop = FALSE], y[ok], se[ok]^2 + 1e-9), y = y)
    }
    main <- model(m)
    feasible <- vapply(observed$.config, function(cfg) meets_constraints(summ[summ$.config == cfg, ], constraints,
                                                                         study$metrics), logical(1))
    incumbent <- if (any(feasible)) max(main$y[feasible], na.rm = TRUE) else min(main$y, na.rm = TRUE)
    fresh <- pool[!key(pool) %in% key(configs), , drop = FALSE]
    if (!nrow(fresh)) break
    xs <- encode_values(study$params, fresh)
    pred <- main$fit(xs)
    score <- expected_improvement(pred$mean, pred$sd, incumbent)
    for (name in names(constraints)) {
      cm <- study$metrics[[name]]
      cp <- model(cm)$fit(xs)
      rule <- constraint_rule(constraints[[name]], cm)
      bound <- metric_sign(cm) * rule$bound
      upper_is_good <- (rule$op %in% c(">=", ">")) == (metric_sign(cm) > 0)
      prob <- if (upper_is_good) stats::pnorm((cp$mean - bound) / cp$sd) else stats::pnorm((bound - cp$mean) / cp$sd)
      score <- if (any(feasible)) score * prob else prob
    }
    order_idx <- order(score, decreasing = TRUE)
    chosen <- NULL
    for (idx in utils::head(order_idx, 50)) {
      values <- lapply(as.list(fresh[idx, ids, drop = FALSE]), function(v) if (is.factor(v)) as.character(v) else v)
      if (is.null(sim_check(object, resolve_args(object, values)))) {
        chosen <- idx
        break
      }
      pool <- pool[key(pool) != key(fresh[idx, , drop = FALSE]), , drop = FALSE]
    }
    if (is.null(chosen)) break
    width <- max(3, nchar(nrow(configs) + iter))
    new_cfg <- tibble::add_column(fresh[chosen, ids, drop = FALSE],
                                  .config = sprintf(paste0("Config%0", width, "d"), nrow(configs) + 1), .before = 1)
    configs <- rbind(configs, new_cfg)
    rows <- c(rows, run_designs(runner, new_cfg, scen, study$resamples))
    raw <- merge_configs(bind_rows_list(rows), configs, ids)
    summ <- summarize_raw(raw, study$metrics, by_scenario = FALSE)
    best <- best_configs(summ, m, constraints, study$metrics, n = 1)
    best_id <- if (nrow(best)) best$.config[1] else NA_character_
    stale <- if (identical(best_id, best_so_far)) stale + 1 else 0
    best_so_far <- best_id
    trace[[length(trace) + 1]] <- tibble::tibble(iteration = it, .config = new_cfg$.config,
                                                 expected_improvement = score[chosen], best = best_id)
    if (isTRUE(control$verbose)) cli::cli_inform("Iteration {it}: tried {new_cfg$.config}; best is {best_id}.")
    if (stale >= no_improve) break
  }
  new_sim_results(object, study, configs, rows, runner, "bayes", if (length(trace)) do.call(rbind, trace) else NULL)
}
