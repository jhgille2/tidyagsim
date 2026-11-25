# Running designs on replicates ------------------------------------------------------------

#' Control the running of a simulation study
#'
#' Like `tune::control_grid()`.
#'
#' @param verbose Print progress.
#' @param parallel Run design batches in parallel with the future framework
#'   (`future.apply` must be installed and a plan set, e.g.
#'   `future::plan(future::multisession, workers = 4)`; this works on Windows).
#' @param cache A directory where each finished replicate is saved. Rerunning
#'   the same study skips them, so an interrupted run resumes where it stopped.
#' @param stop_on_error Stop at the first failing simulation instead of
#'   recording a note and carrying on.
#' @param save_outcomes Keep every outcome the simulator reports (not only the
#'   metrics) in a list column `.outcome`.
#' @return A `control_sim` object.
#' @examples
#' control_sim(verbose = TRUE, cache = tempfile("study-cache"))
#' @export
control_sim <- function(verbose = FALSE, parallel = FALSE, cache = NULL, stop_on_error = FALSE,
                        save_outcomes = FALSE) {
  structure(list(verbose = verbose, parallel = parallel, cache = cache, stop_on_error = stop_on_error,
                 save_outcomes = save_outcomes), class = "control_sim")
}

new_runner <- function(spec, metrics, control) {
  env <- new.env(parent = emptyenv())
  env$spec <- spec
  env$metrics <- metrics
  env$control <- control %||% control_sim()
  env$fingerprint <- spec_fingerprint(spec)
  env$notes <- list()
  env$infeasible <- list()
  if (!is.null(env$control$cache)) dir.create(env$control$cache, showWarnings = FALSE, recursive = TRUE)
  env
}

cache_file <- function(runner, args, seed) {
  if (is.null(runner$control$cache)) return(NULL)
  key <- rlang::hash(list(runner$fingerprint, lapply(args, function(a) if (is.function(a)) deparse(a) else a),
                          as.integer(seed)))
  file.path(runner$control$cache, paste0(key, ".rds"))
}

# configs: tibble with .config and one column per parameter id.
run_designs <- function(runner, configs, scenarios, reps) {
  jobs <- list()
  rows <- list()
  ids <- setdiff(names(configs), ".config")
  for (i in seq_len(nrow(configs))) {
    values <- lapply(as.list(configs[i, ids, drop = FALSE]), function(v) if (is.factor(v)) as.character(v) else v)
    for (sc in scenarios) {
      args <- resolve_args(runner$spec, values, sc$values)
      reason <- sim_check(runner$spec, args)
      if (!is.null(reason)) {
        key <- paste(configs$.config[i], sc$id)
        if (is.null(runner$infeasible[[key]])) {
          runner$infeasible[[key]] <- tibble::tibble(.config = configs$.config[i], .scenario = sc$id,
                                                     reason = reason)
        }
        next
      }
      todo <- integer()
      for (k in seq_len(nrow(reps))) {
        file <- cache_file(runner, args, reps$seed[k])
        if (!is.null(file) && file.exists(file)) {
          rows[[length(rows) + 1]] <- make_row(runner, configs$.config[i], sc$id, reps$id[k], reps$seed[k],
                                               readRDS(file))
        } else {
          todo <- c(todo, k)
        }
      }
      if (length(todo)) {
        jobs[[length(jobs) + 1]] <- list(config = configs$.config[i], scenario = sc$id, args = args,
                                         reps = reps[todo, , drop = FALSE])
      }
    }
  }
  if (!length(jobs)) return(rows)
  total <- sum(vapply(jobs, function(j) nrow(j$reps), integer(1)))
  if (isTRUE(runner$control$verbose)) {
    cli::cli_inform("Simulating {total} replicate{?s} in {length(jobs)} design/scenario batch{?es}.")
  }
  spec <- runner$spec
  stop_on_error <- isTRUE(runner$control$stop_on_error)
  work <- function(job) {
    tryCatch(list(ok = TRUE, value = sim_simulate(spec, job$args, job$reps$seed)),
             error = function(e) {
               if (stop_on_error) stop(e)
               list(ok = FALSE, message = conditionMessage(e))
             })
  }
  results <- if (isTRUE(runner$control$parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) abort("Install future.apply for parallel = TRUE.")
    future.apply::future_lapply(jobs, work, future.seed = TRUE)
  } else {
    out <- vector("list", length(jobs))
    for (j in seq_along(jobs)) {
      out[[j]] <- work(jobs[[j]])
      if (isTRUE(runner$control$verbose) && (j %% 10 == 0 || j == length(jobs))) {
        cli::cli_inform("  {j}/{length(jobs)} batches done")
      }
    }
    out
  }
  for (j in seq_along(jobs)) {
    job <- jobs[[j]]
    res <- results[[j]]
    if (!isTRUE(res$ok)) {
      runner$notes[[length(runner$notes) + 1]] <- tibble::tibble(.config = job$config, .scenario = job$scenario,
                                                                 note = res$message)
      next
    }
    if (length(res$value) != nrow(job$reps)) {
      abort(sprintf("%s returned %d outcomes for %d seeds.", class(spec)[1], length(res$value), nrow(job$reps)))
    }
    for (k in seq_len(nrow(job$reps))) {
      outcome <- res$value[[k]]
      file <- cache_file(runner, job$args, job$reps$seed[k])
      if (!is.null(file)) saveRDS(outcome, file)
      rows[[length(rows) + 1]] <- make_row(runner, job$config, job$scenario, job$reps$id[k], job$reps$seed[k],
                                           outcome)
    }
  }
  rows
}

make_row <- function(runner, config, scenario, id, seed, outcome) {
  values <- lapply(runner$metrics, metric_value, outcome = outcome)
  row <- c(list(.config = config, .scenario = scenario, id = id, seed = as.integer(seed)), values)
  if (isTRUE(runner$control$save_outcomes)) row$.outcome <- list(outcome)
  row
}

bind_rows_list <- function(rows) {
  if (!length(rows)) return(NULL)
  cols <- names(rows[[1]])
  out <- lapply(cols, function(col) {
    v <- lapply(rows, `[[`, col)
    if (col == ".outcome") return(v)
    unlist(v, use.names = FALSE)
  })
  names(out) <- cols
  tibble::as_tibble(out)
}
