# breedsim flowchart programs as a simulation spec ------------------------------------------------

breedsim_modules <- function() {
  rlang::check_installed("reticulate")
  if (!reticulate::py_module_available("breedsim")) {
    abort(c("The Python package `breedsim` was not found.",
            i = "Install it from the plant_breeding_simulator_flowchart repository with `pip install <path>`",
            i = "in the Python that reticulate uses (see reticulate::py_config(), or set RETICULATE_PYTHON)."))
  }
  list(experiment = reticulate::import("breedsim.experiment", convert = FALSE),
       model = reticulate::import("breedsim.model", convert = FALSE),
       schedule = reticulate::import("breedsim.schedule", convert = FALSE),
       simulation = reticulate::import("breedsim.simulation", convert = FALSE),
       comparison = reticulate::import("breedsim.comparison", convert = FALSE),
       rates = reticulate::import("breedsim.rates", convert = FALSE))
}

#' A breedsim breeding program as a simulation spec
#'
#' Wraps a flowchart program drawn in the breedsim desktop editor (or written
#' as JSON) so its fields can be tuned with the tidyagsim verbs. Arguments are
#' named by `"stage.field"` paths, e.g. `` `pyt.reps` = tune() ``. Fields
#' written as `{"tune": ...}` in the program JSON are marked for tuning
#' automatically, with their bounds. Each replicate is simulated by the
#' AlphaSimR engine that ships with breedsim.
#'
#' Reported outcomes per replicate: `release_gain` in `year`; `gain_per_year`
#' and `gain_per_cost` over `start_year`..`year` when `start_year` is given;
#' `expected_heterozygosity` of the crossing pool at `parent_stage` when given;
#' and the program's scheduled `plots`, `genotyped` and `cost` summed over the
#' planning horizon.
#'
#' @param path A breedsim program or study JSON file.
#' @param ... `"stage.field" = value` or `tune()`.
#' @param year Release year scored (defaults to the program horizon).
#' @param start_year First year of the window for rate metrics.
#' @param parent_stage A crossing stage whose parent-pool diversity is reported.
#' @param discount_rate Annual discount rate for `gain_per_cost`.
#' @param settings Path to a breedsim simulation settings JSON file.
#' @param rscript Path to `Rscript` for the engine (defaults to the one on `PATH`
#'   or the `BREEDSIM_RSCRIPT` environment variable).
#' @return A `breeding_program` simulation spec.
#' @examples
#' \dontrun{
#' spec <- breeding_program("examples/conventional.json",
#'                          `pyt.reps` = tune(), `pyt.locations` = tune())
#' res <- tune_grid(spec, sim_replicates(5), grid = 8,
#'                  metrics = sim_metric_set(release_gain, plots))
#' }
#' @export
breeding_program <- function(path, ..., year = NULL, start_year = NULL, parent_stage = NULL, discount_rate = 0,
                             settings = NULL, rscript = NULL) {
  py <- breedsim_modules()
  path <- normalizePath(path, mustWork = TRUE)
  study <- py$experiment$load_study(path)
  program <- reticulate::py_to_r(study$program$model_dump(mode = "json"))
  decisions <- reticulate::py_to_r(study$model_dump(mode = "json"))$decisions
  args <- list()
  bounds <- list()
  for (d in decisions) {
    id <- d$id %||% d$path
    marker <- if (identical(id, d$path)) call("tune") else call("tune", id)
    args[[d$path]] <- marker
    for (link in d$links) args[[link]] <- call("tune", id)
    bounds[[d$path]] <- d[c("low", "high", "choices")]
  }
  spec <- new_sim_spec("breeding_program", args, engine = "alphasimr", program = program, path = path,
                       bounds = bounds, year = year %||% program$years, start_year = start_year,
                       parent_stage = parent_stage, discount_rate = discount_rate,
                       settings = if (is.null(settings)) NULL else normalizePath(settings, mustWork = TRUE),
                       rscript = rscript)
  extra <- list(...)
  if (length(extra)) spec <- set_args.sim_spec(spec, ...)
  spec
}

#' @export
sim_arg_names.breeding_program <- function(object) {
  py <- breedsim_modules()
  prog <- py$model$Program$model_validate(object$program)
  names(reticulate::py_to_r(py$experiment$tunable_fields(prog)))
}

#' @export
sim_parameters.breeding_program <- function(object, ...) {
  py <- breedsim_modules()
  prog <- py$model$Program$model_validate(object$program)
  marked <- tune_args.sim_spec(object)
  out <- list()
  for (path in marked$name) {
    b <- object$bounds[[path]] %||% list()
    decision <- py$experiment$Decision(path = path, low = b$low, high = b$high,
                                       choices = if (is.null(b$choices)) NULL else reticulate::tuple(b$choices))
    exp <- py$experiment$Experiment(program = prog, decisions = reticulate::tuple(decision))
    r <- reticulate::py_get_item(exp$show_ranges(), path)
    get <- function(name) reticulate::py_to_r(reticulate::py_get_attr(r, name))
    kind <- get("kind")
    out[[path]] <- if (kind == "choice") {
      dials::new_qual_param(type = "character", values = unlist(get("choices")), label = stats::setNames(path, path),
                            finalize = NULL)
    } else {
      sim_param(c(get("low"), get("high")), type = if (kind == "integer") "integer" else "double",
                label = sprintf("%s (%s)", path, get("units")))
    }
  }
  out
}

program_with <- function(object, args) {
  prog <- object$program
  ids <- vapply(prog$nodes, function(n) n$id, character(1))
  for (path in names(args)) {
    parts <- strsplit(path, ".", fixed = TRUE)[[1]]
    i <- match(parts[1], ids)
    value <- args[[path]]
    if (is.numeric(value) && length(value) == 1 && value == round(value) && !grepl("cost", parts[2])) {
      value <- as.integer(value)
    }
    if (is.numeric(value) && grepl("cost", parts[2])) value <- format(round(value, 2), nsmall = 2)
    prog$nodes[[i]][[parts[2]]] <- value
  }
  py <- breedsim_modules()
  py$model$Program$model_validate(prog)
}

#' @export
sim_check.breeding_program <- function(object, args, ...) {
  py <- breedsim_modules()
  tryCatch({
    py$schedule$schedule(program_with(object, args))
    NULL
  }, error = function(e) {
    msg <- conditionMessage(e)
    sub("^.*(ScheduleError|ValueError|ValidationError): ?", "", strsplit(msg, "\n")[[1]][1])
  })
}

#' @export
sim_default_metrics.breeding_program <- function(object) {
  if (!is.null(object$start_year)) return(sim_metric_set("gain_per_year", "release_gain", "cost"))
  sim_metric_set("release_gain", "plots", "cost")
}

#' @export
sim_simulate.breeding_program <- function(object, args, seeds, ...) {
  py <- breedsim_modules()
  program <- program_with(object, args)
  settings <- if (is.null(object$settings)) NULL else
    py$simulation$SimulationSettings$model_validate_json(paste(readLines(object$settings, warn = FALSE), collapse = "\n"))
  seeds <- as.list(as.integer(seeds))
  rows <- py$simulation$simulate(program, seeds, settings, object$rscript)
  gains <- reticulate::py_to_r(py$comparison$release_gains(rows, seeds, as.integer(object$year)))
  totals_fn <- reticulate::py_eval(
    "lambda a: dict(plots=sum(r.plots for r in a.values()), genotyped=sum(r.genotyped for r in a.values()), cost=float(sum(r.cost for r in a.values())))",
    convert = TRUE)
  totals <- totals_fn(py$schedule$schedule(program)$annual)
  rates <- NULL
  if (!is.null(object$start_year)) {
    rates <- reticulate::py_to_r(py$rates$gain_rates(program, rows, start_year = as.integer(object$start_year),
                                                     end_year = as.integer(object$year),
                                                     discount_rate = object$discount_rate))
  }
  records <- reticulate::py_to_r(rows)
  lapply(seq_along(seeds), function(k) {
    seed <- seeds[[k]]
    out <- c(list(release_gain = gains[[k]]), totals)
    if (!is.null(rates)) {
      r <- Filter(function(x) x$seed == seed, rates)[[1]]
      out$gain_per_year <- r$gain_per_year
      out$gain_per_cost <- r$gain_per_cost %||% NA_real_
    }
    if (!is.null(object$parent_stage)) {
      pool <- Filter(function(x) x$seed == seed && x$kind == "cross" && identical(x$stage, object$parent_stage) &&
                       x$year == object$year, records)
      out$expected_heterozygosity <- if (length(pool)) mean(vapply(pool, function(x) x$expected_heterozygosity, 0))
                                     else NA_real_
    }
    out
  })
}

#' @export
print.breeding_program <- function(x, ...) {
  cat(sprintf("breeding_program simulation spec (engine: %s)\n  Program: %s (%s years), scoring year %s\n",
              x$engine, x$program$name %||% basename(x$path), x$program$years, x$year))
  marked <- vapply(x$args, is_tune, logical(1))
  if (any(marked)) {
    cat("\nTuned fields:\n")
    for (n in names(x$args)[marked]) cat(sprintf("  %s = %s\n", n, deparse(x$args[[n]])))
  }
  if (any(!marked)) {
    cat("\nChanged fields:\n")
    for (n in names(x$args)[!marked]) cat(sprintf("  %s = %s\n", n, format_arg(x$args[[n]])))
  }
  invisible(x)
}
