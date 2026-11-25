# Precomputes the results shown in vignette("mapping-design").
# Run from the package root: Rscript data-raw/mapping-design.R
# Takes about an hour on one core; set PARALLEL=1 to use future.apply with all cores.
devtools::load_all(".")
if (nzchar(Sys.getenv("PARALLEL"))) {
  future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
}
ctrl <- control_sim(verbose = TRUE, parallel = nzchar(Sys.getenv("PARALLEL")), cache = "data-raw/cache")
out <- list()
t0 <- Sys.time()

# 1. One design --------------------------------------------------------------------------------
ril <- mapping_population(type = "RIL", lines = 200, generation = 6, qtl_pve = 0.05)
reps <- sim_replicates(50)
out$single <- fit_resamples(ril, reps, control = ctrl)

# 2. Tune population size and locations under three QTL effect sizes ---------------------------
ril_grid <- mapping_population(type = "RIL", lines = tune(), locations = tune(), generation = 6)
grid <- expand.grid(lines = c(100L, 200L, 300L, 400L, 500L), locations = 1:3)
scen <- sim_replicates(40, scenarios = sim_scenarios(qtl_pve = c(0.03, 0.05, 0.08)))
out$grid <- tune_grid(ril_grid, scen, grid = grid,
                      metrics = sim_metric_set(power, location_error, false_positive, plots, cost), control = ctrl)

# 3. Compare population types at the same number of lines --------------------------------------
common <- list(lines = 300L, qtl_pve = 0.05)
types <- sim_set(
  F2 = do.call(mapping_population, c(list(type = "F2"), common)),
  RIL = do.call(mapping_population, c(list(type = "RIL", generation = tune()), common)),
  BCsFt = do.call(mapping_population, c(list(type = "BCsFt", backcrosses = tune(), generation = tune()), common)),
  NAM = do.call(mapping_population, c(list(type = "NAM", families = tune(), generation = 6L), common)),
  `AB-NAM` = do.call(mapping_population, c(list(type = "AB-NAM", families = tune(), backcrosses = tune(),
                                                 generation = 6L), common)),
  MAGIC = do.call(mapping_population, c(list(type = "MAGIC", founders = tune(), generation = tune()), common))
)
grids <- list(
  RIL = data.frame(generation = c(3L, 4L, 6L, 8L)),
  BCsFt = expand.grid(backcrosses = 1:2, generation = c(2L, 4L, 6L)),
  NAM = data.frame(families = c(3L, 5L, 10L, 15L, 25L)),
  `AB-NAM` = expand.grid(families = c(5L, 10L, 15L), backcrosses = 1:2),
  MAGIC = expand.grid(founders = c("4", "8"), generation = c(4L, 6L), stringsAsFactors = FALSE)
)
out$types <- sim_map(types, "tune_grid", resamples = sim_replicates(40), grids = grids,
                     metrics = sim_metric_set(power, location_error, support_interval, false_positive, years, cost),
                     control = ctrl, verbose = TRUE)

# 4. Racing a larger space ---------------------------------------------------------------------
ril_race <- mapping_population(type = "RIL", lines = tune(), generation = tune(), qtl_pve = 0.05)
race_params <- extract_parameter_set_dials(ril_race) |>
  update(lines = n_lines(c(100L, 600L)), generation = filial_generation(c(3L, 8L)))
out$race <- tune_race(ril_race, sim_replicates(40), grid = 20, param_info = race_params,
                      metrics = sim_metric_set(power, cost), control = ctrl, burn_in = 5)

# 5. Cheapest design with at least 80% power ---------------------------------------------------
ril_cost <- mapping_population(type = "RIL", lines = tune(), locations = tune(), reps = tune(), generation = 6,
                               qtl_pve = 0.05)
cost_params <- extract_parameter_set_dials(ril_cost) |>
  update(lines = n_lines(c(100L, 800L)), locations = n_locations(c(1L, 4L)), reps = n_reps(c(1L, 2L)))
out$bayes <- tune_bayes(ril_cost, sim_replicates(25), param_info = cost_params,
                        metrics = sim_metric_set(cost, power, plots), constraints = list(power = 0.8),
                        initial = 8, iter = 15, control = ctrl)

out$elapsed <- difftime(Sys.time(), t0, units = "mins")
out$session <- utils::sessionInfo()$R.version$version.string
# Results keep the specs but drop cached outcome lists so the file stays small.
saveRDS(out, "inst/extdata/mapping-design.rds", compress = "xz")
message("Done in ", format(out$elapsed))
