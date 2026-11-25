#' The job-array shard of this process
#'
#' Reads the Slurm array task (`SLURM_ARRAY_TASK_ID`) and the number of tasks
#' (`TIDYAGSIM_TASKS`, set by [write_slurm_array()], or `SLURM_ARRAY_TASK_COUNT`).
#' Outside an array it returns `NULL`, so the same script both runs a shard
#' and, run once more without an array, collects the finished study.
#'
#' @return `c(task, n_tasks)` or `NULL`.
#' @export
sim_shard <- function() {
  task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  n <- Sys.getenv("TIDYAGSIM_TASKS", Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
  if (!nzchar(task) || !nzchar(n)) return(NULL)
  c(as.integer(task), as.integer(n))
}

cluster_profiles <- function() {
  list(
    ceres = list(partition = "ceres", account_required = FALSE,
                 notes = "Every Ceres user has a default account; --account charges a specific project."),
    atlas = list(partition = "atlas", account_required = TRUE,
                 notes = "Atlas rejects jobs without --account (your project name).")
  )
}

#' Write Slurm scripts that split a study across a job array
#'
#' Writes `run.sbatch`, an array job in which each task runs `script` with its
#' shard of the designs, and `collect.sbatch`, which runs after every task
#' succeeds and assembles the results from the shared cache. The study script
#' should build its control with
#' `control_sim(cache = "<shared dir>", shard = sim_shard())` and save the
#' results (e.g. with `saveRDS()`) when it is not a shard
#' (`is.null(sim_shard())`).
#'
#' Profiles exist for SCINet Ceres and Atlas; other clusters work by passing
#' `partition` and `account` directly. Submit with `sbatch run.sbatch`; it
#' submits the collection step itself.
#'
#' @param script Path of the R script to run (copied into `dir`).
#' @param dir Output directory.
#' @param tasks Number of array tasks.
#' @param cluster `"ceres"`, `"atlas"` or `"other"`.
#' @param account Slurm account (required on Atlas).
#' @param time Walltime per task, e.g. `"04:00:00"`.
#' @param memory Memory per task, e.g. `"8G"`.
#' @param cpus CPUs per task (use with `control_sim(parallel = TRUE)`).
#' @param partition Override the profile's partition.
#' @param modules Environment modules to load, e.g. `"r"`; empty when using a container.
#' @param container Optional Apptainer image with R and tidyagsim; tasks then run
#'   `apptainer exec <image> Rscript`.
#' @param max_concurrent Limit on simultaneously running tasks.
#' @return The paths written, invisibly.
#' @examples
#' script <- tempfile(fileext = ".R")
#' writeLines("library(tidyagsim)", script)
#' write_slurm_array(script, tempfile("slurm"), tasks = 20, cluster = "ceres")
#' @export
write_slurm_array <- function(script, dir, tasks, cluster = c("ceres", "atlas", "other"), account = NULL,
                              time = "04:00:00", memory = "4G", cpus = 1, partition = NULL, modules = "r",
                              container = NULL, max_concurrent = NULL) {
  cluster <- match.arg(cluster)
  profile <- cluster_profiles()[[cluster]] %||% list(partition = NULL, account_required = FALSE, notes = NULL)
  if (isTRUE(profile$account_required) && is.null(account)) abort(sprintf("%s needs `account`.", cluster))
  if (!grepl("^(\\d+-)?\\d{1,3}:\\d{2}:\\d{2}$", time)) abort("Use a walltime such as 04:00:00.")
  if (!grepl("^\\d+[KMGT]?$", memory)) abort("Use a memory request such as 8G.")
  if (tasks < 1 || tasks != round(tasks)) abort("`tasks` must be a positive whole number.")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(script, file.path(dir, "study.R"), overwrite = TRUE)
  partition <- partition %||% profile$partition
  header <- function(name, extra = character()) {
    c("#!/bin/bash",
      sprintf("#SBATCH --job-name=%s", name),
      if (!is.null(partition)) sprintf("#SBATCH --partition=%s", partition),
      if (!is.null(account)) sprintf("#SBATCH --account=%s", account),
      sprintf("#SBATCH --time=%s", time),
      sprintf("#SBATCH --mem=%s", memory),
      sprintf("#SBATCH --cpus-per-task=%d", as.integer(cpus)),
      extra,
      "set -euo pipefail",
      "cd \"$SLURM_SUBMIT_DIR\"",
      if (length(modules) && is.null(container)) paste("module load", paste(modules, collapse = " ")),
      if (!is.null(container)) "module load apptainer")
  }
  rscript <- if (is.null(container)) "Rscript" else sprintf("apptainer exec %s Rscript", shQuote(container))
  array <- sprintf("1-%d%s", as.integer(tasks), if (is.null(max_concurrent)) "" else paste0("%", max_concurrent))
  run <- c(header("tidyagsim-run", c(sprintf("#SBATCH --array=%s", array),
                                     "#SBATCH --output=logs/run_%A_%a.out")),
           sprintf("export TIDYAGSIM_TASKS=%d", as.integer(tasks)),
           "if [ \"${SLURM_ARRAY_TASK_ID}\" = \"1\" ]; then",
           "  sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} collect.sbatch",
           "fi",
           sprintf("%s study.R", rscript))
  collect <- c(header("tidyagsim-collect", "#SBATCH --output=logs/collect_%j.out"),
               "unset SLURM_ARRAY_TASK_ID TIDYAGSIM_TASKS",
               sprintf("%s study.R", rscript))
  readme <- c(sprintf("# tidyagsim study for %s", cluster), "",
              "1. Copy this folder to the cluster (e.g. into /project/<your project>/).",
              "2. Submit: `mkdir -p logs && sbatch run.sbatch`",
              sprintf("3. %d array tasks share the designs; the collection job runs when all succeed.", as.integer(tasks)),
              "4. Each finished replicate is saved in the study's cache directory, so resubmitting reruns only what is missing.",
              "", profile$notes)
  paths <- file.path(dir, c("run.sbatch", "collect.sbatch", "README.md"))
  writeLines(run, paths[1])
  writeLines(collect, paths[2])
  writeLines(readme, paths[3])
  invisible(c(file.path(dir, "study.R"), paths))
}

#' Collect every outcome the simulator reported
#'
#' Needs `control_sim(save_outcomes = TRUE)`. Returns one row per design,
#' scenario and replicate, with one column per scalar outcome (for example
#' the QTL position or the number of segregating markers in a mapping
#' population), like `tune::collect_predictions()`.
#'
#' @param x A `sim_results` object.
#' @return A tibble.
#' @export
collect_outcomes <- function(x) {
  check_results(x)
  if (!".outcome" %in% names(x$raw)) abort("Rerun with control_sim(save_outcomes = TRUE).")
  base <- x$raw[, setdiff(names(x$raw), ".outcome"), drop = FALSE]
  outs <- x$raw$.outcome
  keys <- unique(unlist(lapply(outs, function(o) names(Filter(function(v) length(v) == 1 && is.atomic(v), o)))))
  keys <- setdiff(keys, names(base))
  for (k in keys) base[[k]] <- vapply(outs, function(o) if (is.null(o[[k]])) NA_real_ else as.numeric(o[[k]]), 0)
  base
}
