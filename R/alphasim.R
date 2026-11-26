# AlphaSimR engine ------------------------------------------------------------------------------

#' Founder population recipe for AlphaSimR
#'
#' Describes how each replicate builds its founders and trait with
#' [AlphaSimR](https://cran.r-project.org/package=AlphaSimR). Every replicate
#' draws new founders from its own seed, and every design in a study uses the
#' same seeds, so all designs start from the same founders (common random
#' numbers).
#'
#' Any of these fields except `setup` and `haplotypes` can also be given to
#' [alphasim_program()] or [line_pipeline()] as an argument, so they can be
#' tuned or turned into scenarios, for example
#' `sim_scenarios(plot_h2 = c(0.1, 0.3))`.
#'
#' @section Where founder haplotypes come from:
#' * `species = "none"` (the default) draws founder haplotypes with
#'   `AlphaSimR::quickHaplo()`: fast and exactly reproducible from the seed,
#'   but without the linkage disequilibrium a real population history leaves.
#' * `species = "GENERIC"`, `"WHEAT"`, `"MAIZE"`, ... simulates a population
#'   history with MaCS (`AlphaSimR::runMacs()`). MaCS does not follow R's
#'   random seed, so tidyagsim keeps each replicate's founders in memory for
#'   the session: designs in one session share founders, but a new session
#'   (or a parallel worker) draws different ones.
#' * `haplotypes` takes a `MapPop` you made yourself, for example from real
#'   marker data with `AlphaSimR::importHaplo()` or one MaCS run saved with
#'   `saveRDS()`. Each replicate samples `n_founders` of its individuals.
#'
#' @param n_founders Number of inbred founder lines.
#' @param n_chr Number of chromosomes (soybean has 20).
#' @param seg_sites Segregating sites per chromosome.
#' @param qtl_per_chr QTL per chromosome for the additive trait.
#' @param mean,var Founder mean and genetic variance of the trait.
#' @param plot_h2 Heritability of a single plot in the founders; sets the
#'   plot error variance.
#' @param species `"none"` or a MaCS species history; see below.
#' @param haplotypes Optional AlphaSimR `MapPop` to sample founders from.
#' @param setup Optional function of `SP` (the `SimParam`) for anything else,
#'   such as extra traits (`SP$addTraitAD()`), genetic correlations or SNP chips.
#'   When given, it replaces the default trait.
#' @return An `asr_founders` recipe.
#' @examples
#' asr_founders(n_founders = 40, plot_h2 = 0.2)
#' @export
asr_founders <- function(n_founders = 50, n_chr = 20, seg_sites = 200, qtl_per_chr = 20, mean = 0, var = 1,
                         plot_h2 = 0.3, species = "none", haplotypes = NULL, setup = NULL) {
  if (!is.null(haplotypes) && !methods::is(haplotypes, "MapPop")) abort("`haplotypes` must be an AlphaSimR MapPop.")
  structure(list(n_founders = n_founders, n_chr = n_chr, seg_sites = seg_sites, qtl_per_chr = qtl_per_chr,
                 mean = mean, var = var, plot_h2 = plot_h2, species = species, haplotypes = haplotypes,
                 setup = setup),
            class = "asr_founders")
}

#' @export
print.asr_founders <- function(x, ...) {
  if (is.null(x$haplotypes)) {
    src <- if (identical(x$species, "none")) "independent haplotypes" else sprintf("MaCS %s history", x$species)
    cat(sprintf("# AlphaSimR founders: %d inbred lines, %d chromosomes x %d sites (%s)\n",
                as.integer(x$n_founders), as.integer(x$n_chr), as.integer(x$seg_sites), src))
  } else {
    cat(sprintf("# AlphaSimR founders: %d inbred lines sampled from %d supplied haplotypes\n",
                as.integer(x$n_founders), x$haplotypes@nInd))
  }
  cat(sprintf("# Trait: %d QTL per chromosome, plot h2 %s\n", as.integer(x$qtl_per_chr), format(x$plot_h2)))
  invisible(x)
}

founder_fields <- function() setdiff(names(formals(asr_founders)), c("setup", "haplotypes"))

macs_memo <- new.env(parent = emptyenv())

check_alphasimr <- function() {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    abort("This spec needs the AlphaSimR package: install.packages(\"AlphaSimR\").")
  }
}

asr_build <- function(recipe, seed) {
  check_alphasimr()
  fp <- if (!is.null(recipe$haplotypes)) {
    n <- recipe$haplotypes@nInd
    if (recipe$n_founders > n) abort(sprintf("n_founders (%d) is more than the %d supplied haplotypes.",
                                             as.integer(recipe$n_founders), n))
    recipe$haplotypes[sort(sample.int(n, recipe$n_founders))]
  } else if (identical(recipe$species, "none")) {
    AlphaSimR::quickHaplo(nInd = recipe$n_founders, nChr = recipe$n_chr, segSites = recipe$seg_sites,
                          inbred = TRUE)
  } else {
    key <- rlang::hash(list(recipe$n_founders, recipe$n_chr, recipe$seg_sites, recipe$species, seed))
    if (is.null(macs_memo[[key]])) {
      macs_memo[[key]] <- AlphaSimR::runMacs(nInd = recipe$n_founders, nChr = recipe$n_chr,
                                             segSites = recipe$seg_sites, inbred = TRUE, species = recipe$species)
    }
    macs_memo[[key]]
  }
  set.seed(program_seed(as.numeric(seed) + 1))  # trait draws do not depend on how founders were made
  SP <- AlphaSimR::SimParam$new(fp)
  SP$nThreads <- 1L
  if (is.function(recipe$setup)) {
    recipe$setup(SP)
  } else {
    SP$addTraitA(nQtlPerChr = recipe$qtl_per_chr, mean = recipe$mean, var = recipe$var)
    SP$setVarE(h2 = recipe$plot_h2)
  }
  list(pop = AlphaSimR::newPop(fp, simParam = SP), SP = SP)
}

# AlphaSimR functions look up `SP` in the global environment when `simParam`
# is not given. Binding it there for the duration of a replicate lets users
# paste ordinary AlphaSimR scripts; the previous value is restored afterwards.
with_global_sp <- function(SP, expr) {
  genv <- globalenv()
  had <- exists("SP", envir = genv, inherits = FALSE)
  old <- if (had) get("SP", envir = genv)
  assign("SP", SP, envir = genv)
  on.exit(if (had) assign("SP", old, envir = genv) else rm("SP", envir = genv))
  force(expr)
}

#' A simulation spec from an AlphaSimR script
#'
#' Turns an AlphaSimR breeding scheme into a tunable spec. You write one
#' replicate of the scheme as a function of the founders, the `SimParam` and
#' your design arguments; tidyagsim builds the founders for each replicate,
#' seeds the random numbers, runs every design on the same seeds and collects
#' the results. Mark any argument with `tune()` to search over it.
#'
#' Inside `program`, `SP` is also bound globally while it runs, so AlphaSimR
#' calls work with or without `simParam = SP`.
#'
#' Use [asr_tracker()] inside `program` to record the population each year
#' and [asr_summary()] to turn the record into the built-in breeding metrics
#' (`release_gain`, `gain_per_year`, `genetic_variance`, `cost`, ...).
#'
#' @param program A function whose first two arguments are `founders` (an
#'   AlphaSimR `Pop`) and `SP`, followed by design arguments. It returns a
#'   named list of outcomes.
#' @param ... Values for `program`'s arguments and for any [asr_founders()]
#'   field; use `tune()` to mark one for tuning.
#' @param founders An [asr_founders()] recipe.
#' @param params Named list of parameter objects (see [sim_param()]) giving
#'   ranges for tuned arguments.
#' @param check Optional function of the design arguments returning `NULL`
#'   when a design is feasible or a reason string when it is not.
#' @param metrics Optional default [sim_metric_set()].
#' @return An `alphasim_program` spec.
#' @seealso The article "Writing your own AlphaSimR program" builds a program step by step:
#'   \url{https://jhgille2.github.io/tidyagsim/articles/alphasim-programs.html}.
#' @examples
#' \dontrun{
#' recurrent <- alphasim_program(
#'   function(founders, SP, n_select, cycles) {
#'     tr <- asr_tracker(founders)
#'     pop <- founders
#'     for (cycle in seq_len(cycles)) {
#'       pop <- AlphaSimR::setPheno(pop)
#'       parents <- AlphaSimR::selectInd(pop, n_select)
#'       pop <- AlphaSimR::randCross(parents, nCrosses = 100)
#'       tr$record(pop, year = cycle, cost = 100 * 15)
#'     }
#'     asr_summary(tr)
#'   },
#'   n_select = tune(), cycles = 5,
#'   founders = asr_founders(n_chr = 5, seg_sites = 100),
#'   params = list(n_select = sim_param(c(5, 40)))
#' )
#' res <- tune_grid(recurrent, sim_replicates(10), grid = 5)
#' autoplot(res)
#' }
#' @export
alphasim_program <- function(program, ..., founders = asr_founders(), params = list(), check = NULL,
                             metrics = NULL) {
  if (!is.function(program)) abort("`program` must be a function.")
  f <- names(formals(program))
  if (length(f) < 2 || !identical(f[1:2], c("founders", "SP"))) {
    abort("`program` must start with the arguments `founders` and `SP`.")
  }
  if (!inherits(founders, "asr_founders")) abort("`founders` must come from asr_founders().")
  new_sim_spec("alphasim_program", list(...), engine = "AlphaSimR", program = program, founders = founders,
               params = params, check = check, default_metrics = metrics)
}

#' @export
sim_arg_names.alphasim_program <- function(object) {
  f <- names(formals(object$program))
  if ("..." %in% f) return(NULL)
  c(setdiff(f, c("founders", "SP")), founder_fields())
}

#' @export
sim_parameters.alphasim_program <- function(object, ...) {
  c(object$params, founder_parameters())
}

founder_parameters <- function() {
  list(n_founders = sim_param(c(10L, 200L), label = "Founder lines"),
       plot_h2 = plot_heritability(),
       qtl_per_chr = sim_param(c(1L, 200L), label = "QTL per chromosome"))
}

#' @export
sim_check.alphasim_program <- function(object, args, ...) {
  if (is.null(object$check)) return(NULL)
  do.call(object$check, args[intersect(names(args), names(formals(object$check)))])
}

#' @export
sim_default_metrics.alphasim_program <- function(object) {
  object$default_metrics %||% sim_metric_set("release_gain", "gain_per_year", "genetic_variance", "cost")
}

split_founder_args <- function(args, recipe) {
  fields <- intersect(names(args), founder_fields())
  for (f in fields) recipe[[f]] <- args[[f]]
  list(recipe = recipe, args = args[setdiff(names(args), fields)])
}

program_seed <- function(seed) as.integer((as.numeric(seed) * 7919 + 104729) %% 2147483647)

run_alphasim <- function(program, recipe, args, seeds) {
  check_alphasimr()
  lapply(seeds, function(seed) {
    # Founders and the program draw from separate streams, so a MaCS history
    # reused from memory leaves the program's random numbers unchanged.
    built <- with_seed(seed, asr_build(recipe, seed))
    out <- with_seed(program_seed(seed), with_global_sp(
      built$SP, do.call(program, c(list(founders = built$pop, SP = built$SP), args))))
    if (is.data.frame(out)) out <- as.list(out[1, , drop = FALSE])
    as.list(out)
  })
}

#' @export
sim_simulate.alphasim_program <- function(object, args, seeds, ...) {
  s <- split_founder_args(args, object$founders)
  run_alphasim(object$program, s$recipe, s$args, seeds)
}

# Tracking a program over time ---------------------------------------------------------------

#' Record a breeding program as it runs
#'
#' Create a tracker from the founders at the start of an [alphasim_program()],
#' call `$record()` with the population of interest each year (and the money
#' spent that year), then pass the tracker to [asr_summary()]. Genetic values
#' are reported in founder genetic standard deviations, so gains are
#' comparable across traits and scenarios.
#'
#' @param founders The founder `Pop`.
#' @param trait Trait number to follow.
#' @return An `asr_tracker`: a list of functions `record(pop, year, cost = 0, released_lines = NULL)`,
#'   `add_cost(cost)` and `trajectory()`.
#' @export
asr_tracker <- function(founders, trait = 1) {
  check_alphasimr()
  g0 <- AlphaSimR::gv(founders)[, trait]
  base <- list(mean = mean(g0), sd = stats::sd(g0), var = stats::var(g0), he = asr_heterozygosity(founders))
  rows <- list()
  spent <- 0
  released <- NULL
  self <- list(
    base = base,
    record = function(pop, year, cost = 0, released_lines = NULL) {
      spent <<- spent + cost
      g <- AlphaSimR::gv(pop)[, trait]
      rows[[length(rows) + 1]] <<- data.frame(year = year, mean = (mean(g) - base$mean) / base$sd,
                                              genetic_variance = stats::var(g) / base$var,
                                              heterozygosity = asr_heterozygosity(pop) / base$he,
                                              cost = spent)
      if (!is.null(released_lines)) {
        released <<- (AlphaSimR::gv(released_lines)[, trait] - base$mean) / base$sd
      }
      invisible(NULL)
    },
    add_cost = function(cost) {
      spent <<- spent + cost
      invisible(NULL)
    },
    trajectory = function() do.call(rbind, rows),
    released = function() released,
    spent = function() spent
  )
  structure(self, class = "asr_tracker")
}

asr_heterozygosity <- function(pop) {
  geno <- AlphaSimR::pullSegSiteGeno(pop)
  p <- colMeans(geno) / 2
  mean(2 * p * (1 - p))
}

#' Summarize a tracked program into breeding metrics
#'
#' @param tracker An [asr_tracker()].
#' @param ... Extra outcomes to report (for example `plots = 1200`).
#' @return A named list of outcomes: `release_gain` (mean genetic value of
#'   the most recently released lines, or of the last recorded population, in
#'   founder genetic SDs), `gain_per_year`, `genetic_variance` and `expected_heterozygosity`
#'   (last record, relative to the founders), `cost`, `gain_per_cost` (gain per
#'   1000 currency units), `years`, and `trajectory`, a data frame of the
#'   yearly records.
#' @export
asr_summary <- function(tracker, ...) {
  if (!inherits(tracker, "asr_tracker")) abort("`tracker` must come from asr_tracker().")
  traj <- tracker$trajectory()
  if (is.null(traj)) abort("Nothing was recorded; call tracker$record() inside the program.")
  last <- traj[nrow(traj), ]
  gain <- if (length(tracker$released())) mean(tracker$released()) else last$mean
  years <- max(traj$year)
  cost <- tracker$spent()
  c(list(release_gain = gain, gain_per_year = gain / years, genetic_variance = last$genetic_variance,
         expected_heterozygosity = last$heterozygosity, cost = cost,
         gain_per_cost = if (cost > 0) 1000 * gain / cost else NA_real_, years = years),
    list(...), list(trajectory = traj))
}

#' Collect yearly trajectories from a study
#'
#' For specs whose outcomes include a `trajectory` data frame, such as
#' [line_pipeline()] or an [alphasim_program()] that returns [asr_summary()].
#' Needs `control_sim(save_outcomes = TRUE)`.
#'
#' @param x A `sim_results` object.
#' @return A tibble with one row per design, scenario, replicate and year,
#'   including the tuned parameter values.
#' @export
collect_trajectories <- function(x) {
  check_results(x)
  if (!".outcome" %in% names(x$raw)) abort("Rerun with control_sim(save_outcomes = TRUE).")
  keep <- intersect(c(".config", ".scenario", "id", x$params$id), names(x$raw))
  pieces <- lapply(seq_len(nrow(x$raw)), function(i) {
    tr <- x$raw$.outcome[[i]]$trajectory
    if (is.null(tr)) return(NULL)
    info <- as.data.frame(x$raw[rep(i, nrow(tr)), keep, drop = FALSE])
    cbind(info, tr, row.names = NULL)
  })
  out <- do.call(rbind, pieces)
  if (is.null(out)) abort("No outcome has a `trajectory`.")
  tibble::as_tibble(out)
}
