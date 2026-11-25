# A ready-made inbred line pipeline on AlphaSimR ------------------------------------------------

line_defaults <- function() {
  list(n_parents = 30, n_crosses = 50, progeny_per_cross = 20, line_generation = 5, generations_per_year = 2,
       pyt_locations = 1, pyt_reps = 1, n_ayt = 100, ayt_locations = 4, ayt_reps = 2, n_release = 5,
       cycles = 4, gxe_var = 0.5, cross_cost = 25, ssd_cost = 2, plot_cost = 15,
       max_plots = Inf, max_cost = Inf)
}

#' A soybean-style inbred line breeding pipeline
#'
#' A ready-made [alphasim_program()] for the pipeline most public soybean (and
#' other self-pollinated crop) programs run: cross elite parents, advance the
#' progeny to inbred lines by single-seed descent, test every line in a
#' preliminary yield trial (PYT), advance the best to an advanced yield trial
#' (AYT), release the best AYT lines and recycle the best as parents for the
#' next cycle. Any argument can be tuned with `tune()`.
#'
#' @section How one cycle works:
#' 1. `n_crosses` random crosses among the parents (the founders in cycle 1).
#' 2. Each F1 gives `progeny_per_cross` F2 plants, advanced by single-seed
#'    descent to F`line_generation` lines, so a cycle makes
#'    `n_crosses * progeny_per_cross` lines.
#' 3. PYT: every line at `pyt_locations` x `pyt_reps` plots. The top `n_ayt`
#'    lines by trial mean go on.
#' 4. AYT: `ayt_locations` x `ayt_reps` plots. The top `n_release` lines are
#'    released and the top `n_parents` become the next cycle's parents.
#'
#' Trial means carry error variance
#' `gxe_var / locations + plot_var / (locations * reps)`, where `plot_var`
#' comes from the founders' `plot_h2` and `gxe_var` is in units of the founder
#' genetic variance, so extra locations average out genotype-by-environment
#' noise that extra reps cannot.
#'
#' @section Time and cost:
#' A cycle takes `line_generation / generations_per_year` years in the
#' nursery (crossing and selfing) plus one year for each trial. Cycles do not
#' overlap, which keeps the model simple but understates the gain per year of
#' programs that start new crosses every year. Cost adds crosses
#' (`cross_cost`), single-seed-descent plants per generation (`ssd_cost`) and
#' trial plots (`plot_cost`). Designs whose plots or cost per cycle exceed
#' `max_plots` or `max_cost` are skipped as infeasible.
#'
#' @section Outcomes:
#' `release_gain` (genetic value of the last cycle's releases, in founder
#' genetic SDs), `gain_per_year`, `gain_per_cost` (per 1000 currency units),
#' `genetic_variance` and `expected_heterozygosity` of the final parents
#' (relative to the founders), `accuracy` and `ayt_accuracy` (correlation of
#' PYT and AYT means with genetic value, among the lines in each trial), `lines`, `plots` and `cost_per_cycle` per cycle, total
#' `cost`, `years`, and a
#' yearly `trajectory` (see [collect_trajectories()]).
#'
#' @param ... Arguments to set or tune: `n_parents`, `n_crosses`,
#'   `progeny_per_cross`, `line_generation`, `generations_per_year`,
#'   `pyt_locations`, `pyt_reps`, `n_ayt`, `ayt_locations`, `ayt_reps`,
#'   `n_release`, `cycles`, `gxe_var`, `cross_cost`, `ssd_cost`, `plot_cost`,
#'   `max_plots`, `max_cost`, and any [asr_founders()] field such as `plot_h2`.
#' @param founders An [asr_founders()] recipe.
#' @return A `line_pipeline` spec.
#' @examples
#' line_pipeline(n_ayt = tune(), ayt_locations = tune())
#' \dontrun{
#' spec <- line_pipeline(n_ayt = tune(), ayt_locations = tune(), cycles = 3)
#' res <- tune_grid(spec, sim_replicates(10), grid = 9)
#' show_best(res, metric = "gain_per_year")
#' }
#' @export
line_pipeline <- function(..., founders = asr_founders()) {
  args <- list(...)
  if (!inherits(founders, "asr_founders")) abort("`founders` must come from asr_founders().")
  defaults <- line_defaults()
  defaults[names(args)] <- args
  spec <- new_sim_spec("line_pipeline", defaults, engine = "AlphaSimR", founders = founders)
  check_arg_names(spec, names(args))
  spec
}

#' @export
sim_arg_names.line_pipeline <- function(object) c(names(line_defaults()), founder_fields())

#' @export
sim_parameters.line_pipeline <- function(object, ...) {
  c(list(n_parents = sim_param(c(5L, 100L), label = "Parents recycled"),
         n_crosses = sim_param(c(10L, 200L), label = "Crosses per cycle"),
         progeny_per_cross = sim_param(c(5L, 100L), label = "Lines per cross"),
         line_generation = filial_generation(c(3L, 8L)),
         generations_per_year = generations_per_year(),
         pyt_locations = n_locations(), pyt_reps = n_reps(),
         n_ayt = sim_param(c(10L, 300L), label = "Lines advanced to AYT"),
         ayt_locations = n_locations(c(1L, 12L)), ayt_reps = n_reps(),
         n_release = sim_param(c(1L, 20L), label = "Lines released per cycle"),
         cycles = sim_param(c(1L, 10L), label = "Breeding cycles"),
         gxe_var = sim_param(c(0, 2), type = "double", label = "GxE variance (founder Vg units)"),
         plot_cost = sim_param(c(5, 50), type = "double", label = "Cost per plot")),
    founder_parameters())
}

#' @export
sim_default_metrics.line_pipeline <- function(object) {
  sim_metric_set("gain_per_year", "release_gain", "genetic_variance", "accuracy", "plots", "cost")
}

line_resources <- function(a) {
  lines <- a$n_crosses * a$progeny_per_cross
  plots <- lines * a$pyt_locations * a$pyt_reps + a$n_ayt * a$ayt_locations * a$ayt_reps
  ssd_plants <- lines * max(a$line_generation - 2, 0)
  cost <- a$n_crosses * a$cross_cost + ssd_plants * a$ssd_cost + plots * a$plot_cost
  years <- a$line_generation / a$generations_per_year + 2
  list(lines = lines, plots = plots, cost = cost, years = years)
}

#' @export
sim_check.line_pipeline <- function(object, args, ...) {
  a <- args
  r <- line_resources(a)
  if (a$n_ayt > r$lines) return(sprintf("n_ayt (%d) is more than the %d lines made.", a$n_ayt, r$lines))
  if (a$n_parents > a$n_ayt) return("n_parents is more than n_ayt.")
  if (a$n_release > a$n_ayt) return("n_release is more than n_ayt.")
  if (a$n_parents < 2) return("Use at least 2 parents.")
  if (r$plots > a$max_plots) return(sprintf("%d plots per cycle is over max_plots.", r$plots))
  if (r$cost > a$max_cost) return(sprintf("Cost per cycle %.0f is over max_cost.", r$cost))
  NULL
}

line_program <- function(founders, SP, n_parents, n_crosses, progeny_per_cross, line_generation,
                         generations_per_year, pyt_locations, pyt_reps, n_ayt, ayt_locations, ayt_reps,
                         n_release, cycles, gxe_var, cross_cost, ssd_cost, plot_cost, max_plots, max_cost) {
  a <- as.list(environment())
  r <- line_resources(a)
  tr <- asr_tracker(founders)
  plot_var <- as.numeric(SP$varE)[1]
  gxe <- gxe_var * tr$base$var
  trial_var <- function(locs, reps) gxe / locs + plot_var / (locs * reps)
  parents <- founders
  year <- 0
  acc <- numeric()
  ayt_acc <- numeric()
  for (cycle in seq_len(cycles)) {
    f1 <- AlphaSimR::randCross(parents, nCrosses = n_crosses, nProgeny = 1, simParam = SP)
    lines <- AlphaSimR::self(f1, nProgeny = progeny_per_cross, simParam = SP)
    for (g in seq_len(max(line_generation - 2, 0))) lines <- AlphaSimR::self(lines, simParam = SP)
    pyt <- AlphaSimR::setPheno(lines, varE = trial_var(pyt_locations, pyt_reps), simParam = SP)
    acc <- c(acc, suppressWarnings(stats::cor(AlphaSimR::gv(pyt)[, 1], AlphaSimR::pheno(pyt)[, 1])))
    ayt <- AlphaSimR::selectInd(pyt, n_ayt, use = "pheno", simParam = SP)
    ayt <- AlphaSimR::setPheno(ayt, varE = trial_var(ayt_locations, ayt_reps), simParam = SP)
    ayt_acc <- c(ayt_acc, suppressWarnings(stats::cor(AlphaSimR::gv(ayt)[, 1], AlphaSimR::pheno(ayt)[, 1])))
    released <- AlphaSimR::selectInd(ayt, n_release, use = "pheno", simParam = SP)
    parents <- AlphaSimR::selectInd(ayt, n_parents, use = "pheno", simParam = SP)
    year <- year + r$years
    tr$record(parents, year = year, cost = r$cost, released_lines = released)
  }
  asr_summary(tr, accuracy = mean(acc, na.rm = TRUE), ayt_accuracy = mean(ayt_acc, na.rm = TRUE), lines = r$lines, plots = r$plots,
              crosses = n_crosses, cost_per_cycle = r$cost)
}

#' @export
sim_simulate.line_pipeline <- function(object, args, seeds, ...) {
  s <- split_founder_args(args, object$founders)
  run_alphasim(line_program, s$recipe, s$args, seeds)
}
